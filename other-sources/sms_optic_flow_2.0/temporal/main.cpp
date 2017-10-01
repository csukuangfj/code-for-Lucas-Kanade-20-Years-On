// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#include "brox_optic_flow.h"

#include <iostream>
#include <cstdio>

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
#define PAR_DEFAULT_DIR "./"
#define PAR_DEFAULT_VERBOSE 0


using namespace std;


/**
 *
 *  Function to read images using the iio library
 *  It allocates memory for the image and returns true if it
 *  correctly reads the image.
 *
 */
bool read_image(const char *fname, float **f, int *nx, int *ny)
{
    *f = iio_read_image_float(fname, nx, ny);

    return *f ? true : false;
}


/**
 *
 *  Function to read all the images of the sequence
 *
 */
bool read_images(char *argv[], int frames, float **I, int &nx, int &ny)
{
    bool correct = true;

    int pos = 0;
    int i   = 0;

    //read all the images in the sequence
    while(correct && i < frames)
    {
	float *f;
	int nxx, nyy;
	
	//read the next image
	correct = read_image(argv[i], &f, &nxx, &nyy);
	
	//if it is the first image, allocate memory
	if(correct && i == 0)
	{
	    nx = nxx; ny = nyy;
	    *I = new float[nx * ny * frames];
	}
	else
	{
	    correct = correct && nx == nxx && ny == nyy;
	}
	
	//if the read image is correct, store it in the array
	if(correct)
	    for(int j = 0; j < nx * ny; j++)
		(*I)[pos++] = f[j];

	free(f);
	
	i++;
    }

    return correct;
}


/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -nimages     number of images in the sequence
 *   -I1          first image
 *   -...         sequence of images
 *   -In          last image
 *   -alpha       smoothing parameter
 *   -gamma       gradient constancy parameter
 *   -nscales     number of scales for the pyramidal approach
 *   -zoom_factor reduction factor for creating the scales
 *   -TOL         stopping criterion threshold for the iterative process
 *   -inner_iter  number of inner iterations
 *   -outer_iter  number of outer iterations
 *   -dir         output directory to store the flows (flow00.uv ... flowNN.uv)
 *   -verbose     switch on/off messages
 *
 */
int main(int argc, char *argv[])
{
    if(argc < 3)

	cout << "Usage: " << argv[0]
	      << " nimages I1...In"
	      << " [alpha gamma nscales zoom_factor"
	      <<	 " TOL inner_iter outer_iter dir verbose]"
	      << endl;
	      
    else
    {
	int i = 1;
	int nx=-1, ny=-1;
	float *I=0; 

	//read the parameters
	const int frames  = atoi(argv[i]); i++;

	//read the input images
	bool correct  = read_images(&argv[i], frames, &I, nx, ny); i += frames;

	//read more parameters
	float alpha   = (argc > i)? atof(argv[i]):PAR_DEFAULT_ALPHA; i++;
	float gamma   = (argc > i)? atof(argv[i]):PAR_DEFAULT_GAMMA; i++;
	int   nscales = (argc > i)? atoi(argv[i]):PAR_DEFAULT_NSCALES; i++;
	float zfactor = (argc > i)? atof(argv[i]):PAR_DEFAULT_ZFACTOR; i++;
	float TOL     = (argc > i)? atof(argv[i]):PAR_DEFAULT_TOL; i++;
	int   initer  = (argc > i)? atoi(argv[i]):PAR_DEFAULT_INNER_ITER; i++;
	int   outiter = (argc > i)? atoi(argv[i]):PAR_DEFAULT_OUTER_ITER; i++;
	const char *dir = (argc > i)? argv[i]:PAR_DEFAULT_DIR; i++;
	int   verbose = (argc > i)? atoi(argv[i]):PAR_DEFAULT_VERBOSE; i++;

	//check parameters
	if(alpha   <= 0) alpha   = PAR_DEFAULT_ALPHA;
	if(gamma   <  0) gamma   = PAR_DEFAULT_GAMMA;
	if(nscales <= 0) nscales = PAR_DEFAULT_NSCALES;
	if(zfactor <= 0 || zfactor >= 1) zfactor = PAR_DEFAULT_ZFACTOR;
	if(TOL     <= 0) TOL = PAR_DEFAULT_TOL;
	if(initer  <= 0) initer  = PAR_DEFAULT_INNER_ITER;
	if(outiter <= 0) outiter = PAR_DEFAULT_OUTER_ITER;

	// if the images are correct, compute the optical flow
	if(correct)
	{  
	    //set the number of scales according to the size of the
	    //images.  The value N is computed to assure that the smaller
	    //images of the pyramid don't have a size smaller than 16x16
	    const float N = 1 + log(std::min(nx, ny) / 16.) / log(1./zfactor);
	    if((int) N < nscales) nscales = (int) N;

	    if(verbose)   
		cout  << endl <<  " alpha:" << alpha << " gamma:" << gamma 
		      << " scales:" << nscales << " nu:" << zfactor 
		      << " TOL:" << TOL << " inner:" << initer << " outer:" << outiter 
		      << endl;
		  
	    //allocate memory for the flow
	    float *u = new float[nx*ny*(frames-1)];
	    float *v = new float[nx*ny*(frames-1)];

	    //compute the optic flows 
	    brox_optic_flow(
		I, u, v, nx, ny, frames, alpha, gamma, 
		nscales, zfactor, TOL, initer, outiter, verbose
	    );

	    //save the flows in directory
	    for(int i = 0; i < frames-1; i++)
	    {
		char file[30];
		sprintf(file, "%s/flow%.2d.uv", dir, i);
		float *f = new float[nx * ny * 2];
		for (int j = 0; j < nx * ny; j++) {
		    f[2*j]   = u[i*nx*ny+j];
		    f[2*j+1] = v[i*nx*ny+j];
		}
		iio_save_image_float_vec(file, f, nx, ny, 2);
		delete []f;
	    }

	    //free dynamic memory
	    delete []I;
	    delete []u;
	    delete []v;
	}
	else
	    cerr << "Cannot read the images or the size of the images are not equal" << endl;
    }

    return 0;
}

