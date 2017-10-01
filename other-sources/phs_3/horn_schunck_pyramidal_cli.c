// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef DISABLE_OMP
#include <omp.h>
#endif//DISABLE_OMP


#include "iio.h"

#include "xmalloc.c"
#include "horn_schunck_pyramidal.c"



#define PAR_DEFAULT_NPROC 0
#define PAR_DEFAULT_ALPHA 7
#define PAR_DEFAULT_NSCALES 10
#define PAR_DEFAULT_ZFACTOR 0.5
#define PAR_DEFAULT_NWARPS 10
#define PAR_DEFAULT_TOL 0.0001
#define PAR_DEFAULT_MAXITER 150
#define PAR_DEFAULT_VERBOSE 0
#define PAR_MAX_ZFACTOR 0.99



/**
 *
 *  Function to read images using the iio library
 *  It allocates memory for the image and returns true if it
 *  correctly reads the image.
 *
 */
static int read_image(
	const char *fname, // file name of the image
	float **f,         // allocated memory for the image
	int *w,            // image width
	int *h             // image height
)
{
	*f = iio_read_image_float(fname, w, h);
	return *f ? true : false;
}



/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -I1          first image
 *   -I2          second image
 *   -out_file    name of the output flow field
 *   -processors  number of threads used with the OpenMP library
 *   -alpha       smoothing parameter
 *   -nscales     number of scales for the pyramidal approach
 *   -zoom_factor reduction factor for creating the scales
 *   -nwarps      number of warps per scales
 *   -TOL         stopping criterion threshold for the iterative process
 *   -maxiter     maximum number of iterations
 *   -verbose     switch on/off messages
 *
 */
int main(int argc, char *argv[])
{
	if(argc < 3) {
		fprintf(stderr, "Usage: %s"
				" I1 I2 [out_file processors"
				" alpha nscales zoom_factor"
				" nwarps TOL maxiter verbose]\n", *argv);
		return EXIT_FAILURE;
	}

	int i = 1;

	//read the parameters and set default values
	const char *image1  = argv[i]; i++;
	const char *image2  = argv[i]; i++;
	char *outfile = (argc >= 4)?  argv[i]:"flow.flo"; i++;

	int   nproc   = (argc >= 5)?  atoi(argv[i]):PAR_DEFAULT_NPROC;   i++;
	float alpha   = (argc >= 6)?  atof(argv[i]):PAR_DEFAULT_ALPHA;   i++;
	int   nscales = (argc >= 7)?  atoi(argv[i]):PAR_DEFAULT_NSCALES; i++;
	float zfactor = (argc >= 8)?  atof(argv[i]):PAR_DEFAULT_ZFACTOR; i++;
	int   warps   = (argc >= 9)?  atoi(argv[i]):PAR_DEFAULT_NWARPS;  i++;
	float TOL     = (argc >= 10)? atof(argv[i]):PAR_DEFAULT_TOL;     i++;
	int   maxiter = (argc >= 11)? atoi(argv[i]):PAR_DEFAULT_MAXITER; i++;
	int   verbose = (argc >= 12)? atoi(argv[i]):PAR_DEFAULT_VERBOSE; i++;

	// check parameters
#ifndef DISABLE_OMP
	if(nproc   >  0) omp_set_num_threads(nproc);
#endif//DISABLE_OMP
	if(alpha   <= 0) alpha   = PAR_DEFAULT_ALPHA;
	if(nscales <= 0) nscales = PAR_DEFAULT_NSCALES;
	if(zfactor <= 0) zfactor = PAR_DEFAULT_ZFACTOR;
	if(zfactor >= 1) zfactor = PAR_MAX_ZFACTOR;
	if(warps   <= 0) warps   = PAR_DEFAULT_NWARPS;
	if(TOL     <= 0) TOL     = PAR_DEFAULT_TOL;

	float *I1, *I2;
	int    nx, ny, nx1, ny1;

	// read the input images
	bool correct1 = read_image(image1, &I1, &nx, &ny);
	bool correct2 = read_image(image2, &I2, &nx1, &ny1);

	// if the images are correct, compute the optical flow
	if(correct1 && correct2 && nx == nx1 && ny == ny1)
	{
		// Set the number of scales according to the size of the
		// images.  The value N is computed to assure that the smaller
		// images of the pyramid don't have a size smaller than 16x16
		const float N = 1 + log(hypot(nx, ny)/16) / log(1/zfactor);
		if(N < nscales)
			nscales = N;

		if (verbose)
			fprintf(stderr, "nproc=%d alpha=%g nscales=%d "
					"zfactor=%g warps=%d epsilon=%g\n",
					nproc, alpha, nscales, zfactor, warps,
					TOL);

		// allocate memory for the optical flow
		float *u = xmalloc(2 * nx * ny * sizeof*u);
		float *v = u + nx*ny;

		// compute the optical flow
		horn_schunck_pyramidal(
				I1, I2, u, v, nx, ny,
				alpha, nscales, zfactor, warps, TOL, maxiter,
				verbose
				);

		// save the flow
		iio_save_image_float_split(outfile, u, nx, ny, 2);

		// free allocated memory
		free(I1);
		free(I2);
		free(u);
	} else {
		fprintf(stderr, "Cannot read the input images "
				"or their sizes are different.\n");
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
