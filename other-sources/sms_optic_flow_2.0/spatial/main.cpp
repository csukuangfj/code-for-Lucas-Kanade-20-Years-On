// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#include <algorithm>
#include <iostream>

#include "brox_optic_flow.h"


extern "C" {
#include "iio.h"
}

#define PAR_DEFAULT_ALPHA 18
#define PAR_DEFAULT_GAMMA 7
#define PAR_DEFAULT_NSCALES 100
#define PAR_DEFAULT_ZFACTOR 0.75
#define PAR_DEFAULT_TOL 0.0001
#define PAR_DEFAULT_INNER_ITER 1
#define PAR_DEFAULT_OUTER_ITER 15
#define PAR_DEFAULT_VERBOSE 0


using namespace std;


/**
 *
 *  Function to read images using the iio library
 *  It allocates memory for the image and returns true if it
 *  correctly reads the image.
 *
 */
bool read_image(const char *fname, float **f, int *w, int *h)
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
 *   -alpha       smoothing parameter
 *   -gamma       gradient constancy parameter
 *   -nscales     number of scales for the pyramidal approach
 *   -zoom_factor reduction factor for creating the scales
 *   -TOL         stopping criterion threshold for the iterative process
 *   -inner_iter  number of inner iterations
 *   -outer_iter  number of outer iterations
 *   -verbose     switch on/off messages
 *
 */
int main(int argc, char *argv[])
{
    if(argc < 3)

	cout << "Usage: " << argv[0]
	      << " I1 I2 [out_file"
	      << " alpha gamma nscales zoom_factor"
	      <<	" TOL inner_iter outer_iter verbose]"
	      << endl;
	      
    else
    {
	int i = 1;

	//read parameters from the console
	const char *image1  = argv[i]; i++;
	const char *image2  = argv[i]; i++;
	const char *outfile = (argc >= 4)?  argv[i]:"flow.uv"; i++;

	float alpha   = (argc > i)? atof(argv[i]):PAR_DEFAULT_ALPHA; i++;
	float gamma   = (argc > i)? atof(argv[i]):PAR_DEFAULT_GAMMA; i++;
	int   nscales = (argc > i)? atoi(argv[i]):PAR_DEFAULT_NSCALES; i++;
	float zfactor = (argc > i)? atof(argv[i]):PAR_DEFAULT_ZFACTOR; i++;
	float TOL     = (argc > i)? atof(argv[i]):PAR_DEFAULT_TOL; i++;
	int   initer  = (argc > i)? atoi(argv[i]):PAR_DEFAULT_INNER_ITER; i++;
	int   outiter = (argc > i)? atoi(argv[i]):PAR_DEFAULT_OUTER_ITER; i++;
	int   verbose = (argc > i)? atoi(argv[i]):PAR_DEFAULT_VERBOSE; i++;

	//check parameters
	if(alpha   <= 0) alpha   = PAR_DEFAULT_ALPHA;
	if(gamma   <  0) gamma   = PAR_DEFAULT_GAMMA;
	if(nscales <= 0) nscales = PAR_DEFAULT_NSCALES;
	if(zfactor <= 0 || zfactor >= 1) zfactor = PAR_DEFAULT_ZFACTOR;
	if(TOL     <= 0) TOL = PAR_DEFAULT_TOL;
	if(initer  <= 0) initer  = PAR_DEFAULT_INNER_ITER;
	if(outiter <= 0) outiter = PAR_DEFAULT_OUTER_ITER;

	int nx, ny, nx1, ny1;
	
	float *I1, *I2; 
	
	//read the input images
	bool correct1 = read_image(image1, &I1, &nx, &ny);
	bool correct2 = read_image(image2, &I2, &nx1, &ny1);

	// if the images are correct, compute the optical flow
	if(correct1 && correct2 && nx == nx1 && ny == ny1)
	{  
	    //set the number of scales according to the size of the
	    //images.  The value N is computed to assure that the smaller
	    //images of the pyramid don't have a size smaller than 16x16
	    const float N = 1 + log(std::min(nx, ny) / 16.) / log(1./zfactor);
	    if((int) N < nscales) nscales = (int) N;

	    if(verbose)
		cout  << endl << " alpha:" << alpha << " gamma:" << gamma
		      << " scales:" << nscales << " nu:" << zfactor 
		      << " TOL:" << TOL << " inner:" << initer << " outer:" << outiter 
		      << endl;
		  
	    //allocate memory for the flow
	    float *u = new float[nx*ny];
	    float *v = new float[nx*ny];

	    //compute the optic flow
	    brox_optic_flow(
		I1, I2, u, v, nx, ny, alpha, gamma, 
		nscales, zfactor, TOL, initer, outiter, verbose
	    );
	    
	    //save the flow 
	    float *f = new float[nx * ny * 2];
	    for (int i = 0; i < nx * ny; i++) {
		f[2*i] = u[i];
		f[2*i+1] = v[i];
	    }
	    iio_save_image_float_vec((char *)outfile, f, nx, ny, 2);

	    //free dynamic memory
	    free(I1);
	    free(I2);
	    delete []u;
	    delete []v;
	    delete []f;
	}
	else
	    cerr << "Cannot read the images or the size of the images are not equal" << endl;
    }

    return 0;
}
