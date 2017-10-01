// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef MASK_C
#define MASK_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "xmalloc.c"


#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING

/**
 *
 * Compute the gradient of an image using centered differences
 *
 */
void gradient(
	const float *input, // input image
	float *dx,          // computed x derivative
	float *dy,          // computed y derivative
	const int nx,       // image width
	const int ny        // image height
)
{
	// compute gradient in the central body of the image
	#pragma omp parallel for
	for(int i = 1; i < ny-1; i++)
	{
		for(int j = 1; j < nx-1; j++)
		{
			const int k = i * nx + j;
			dx[k] = 0.5*(input[k+1] - input[k-1]);
			dy[k] = 0.5*(input[k+nx] - input[k-nx]);
		}
	}

	// compute gradient in the first and last rows
	for(int j = 1; j < nx-1; j++)
	{
		dx[j] = 0.5*(input[j+1] - input[j-1]);
		dy[j] = 0.5*(input[j+nx] - input[j]);

		const int k = (ny - 1) * nx + j;

		dx[k] = 0.5*(input[k+1] - input[k-1]);
		dy[k] = 0.5*(input[k] - input[k-nx]);
	}

	// compute gradient in the first and last columns
	for(int i = 1; i < ny-1; i++)
	{
		const int p = i * nx;
		dx[p] = 0.5*(input[p+1] - input[p]);
		dy[p] = 0.5*(input[p+nx] - input[p-nx]);

		const int k = (i+1) * nx - 1;

		dx[k] = 0.5*(input[k] - input[k-1]);
		dy[k] = 0.5*(input[k+nx] - input[k-nx]);
	}

	// compute the gradient in the corners
	dx[0] = 0.5*(input[1] - input[0]);
	dy[0] = 0.5*(input[nx] - input[0]);

	dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
	dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

	dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
	dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

	dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
	dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);

}


/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
void gaussian(
	float *I,             // input/output image
	const int xdim,       // image width
	const int ydim,       // image height
	const double sigma    // Gaussian sigma
)
{
	const int boundary_condition = DEFAULT_BOUNDARY_CONDITION;
	const int window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

	const double den  = 2*sigma*sigma;
	const int   size = (int) (window_size * sigma) + 1 ;
	const int   bdx  = xdim + size;
	const int   bdy  = ydim + size;

	if (boundary_condition && size > xdim) {
		fprintf(stderr, "GaussianSmooth: sigma too large\n");
		abort();
	}

	// compute the coefficients of the 1D convolution kernel
	double B[size];
	for(int i = 0; i < size; i++)
		B[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) * exp(-i * i / den);

	// normalize the 1D convolution kernel
	double norm = 0;
	for(int i = 0; i < size; i++)
		norm += B[i];
	norm *= 2;
	norm -= B[0];
	for(int i = 0; i < size; i++)
		B[i] /= norm;

	// convolution of each line of the input image
	double *R = xmalloc((size + xdim + size)*sizeof*R);

	for (int k = 0; k < ydim; k++)
	{
		int i, j;
		for (i = size; i < bdx; i++)
			R[i] = I[k * xdim + i - size];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for(i = 0, j = bdx; i < size; i++, j++)
				R[i] = R[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + size-i];
				R[j] = I[k * xdim + xdim-i-1];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for(i = 0, j = bdx; i < size; i++, j++) {
				R[i] = I[k * xdim + xdim-size+i];
				R[j] = I[k * xdim + i];
			}
			break;
		}

		for (i = size; i < bdx; i++)
		{
			double sum = B[0] * R[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * ( R[i-j] + R[i+j] );
			I[k * xdim + i - size] = sum;
		}
	}

	// convolution of each column of the input image
	double *T = xmalloc((size + ydim + size)*sizeof*T);

	for (int k = 0; k < xdim; k++)
	{
		int i, j;
		for (i = size; i < bdy; i++)
			T[i] = I[(i - size) * xdim + k];

		switch (boundary_condition)
		{
		case BOUNDARY_CONDITION_DIRICHLET:
			for (i = 0, j = bdy; i < size; i++, j++)
				T[i] = T[j] = 0;
			break;

		case BOUNDARY_CONDITION_REFLECTING:
			for (i = 0, j = bdy; i < size; i++, j++) {
				T[i] = I[(size-i) * xdim + k];
				T[j] = I[(ydim-i-1) * xdim + k];
			}
			break;

		case BOUNDARY_CONDITION_PERIODIC:
			for( i = 0, j = bdx; i < size; i++, j++) {
				T[i] = I[(ydim-size+i) * xdim + k];
				T[j] = I[i * xdim + k];
			}
			break;
		}

		for (i = size; i < bdy; i++)
		{
			double sum = B[0] * T[i];
			for (j = 1; j < size; j++ )
				sum += B[j] * (T[i-j] + T[i+j]);
			I[(i - size) * xdim + k] = sum;
		}
	}

	free(R);
	free(T);
}


#endif//MASK_C
