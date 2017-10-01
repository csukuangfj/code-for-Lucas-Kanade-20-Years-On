// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef HORN_SCHUNCK_PYRAMIDAL_C
#define HORN_SCHUNCK_PYRAMIDAL_C

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mask.c"
#include "zoom.c"
#include "xmalloc.c"


#define SOR_EXTRAPOLATION_PARAMETER 1.9
#define INPUT_PRESMOOTHING_SIGMA 0.8



/**
 *
 *  Function to compute the SOR iteration at a given position
 *  (SOR = Successive Over-Relaxation)
 *
 */
static float sor_iteration(
	const float *Au, // constant part of the numerator of u
	const float *Av, // constant part of the numerator of v
	const float *Du, // denominator of u
	const float *Dv, // denominator of v
	const float *D,  // constant part of the numerator
	float 	    *u,  // x component of the flow
	float 	    *v,  // y component of the flow
	const float  al, // alpha smoothness parameter
	const int    p,  // current position
	const int    p1, // up-left neighbor
	const int    p2, // up-right neighbor
	const int    p3, // bottom-left neighbor
	const int    p4, // bottom-right neighbor
	const int    p5, // up neighbor
	const int    p6, // left neighbor
	const int    p7, // bottom neighbor
	const int    p8  // right neighbor
)
{
	// set the SOR extrapolation parameter
	const float w = SOR_EXTRAPOLATION_PARAMETER;

	// compute the divergence
	const float ula = 1./12. * (u[p1] + u[p2] + u[p3] + u[p4]) +
			  1./6.  * (u[p5] + u[p6] + u[p7] + u[p8]);
	const float vla = 1./12. * (v[p1] + v[p2] + v[p3] + v[p4]) +
			  1./6.  * (v[p5] + v[p6] + v[p7] + v[p8]);

	// store the previous values
	const float uk = u[p];
	const float vk = v[p];

	// update the flow
	u[p] = (1.0 - w) * uk + w * (Au[p] - D[p] * v[p] + al * ula) / Du[p];
	v[p] = (1.0 - w) * vk + w * (Av[p] - D[p] * u[p] + al * vla) / Dv[p];

	// return the convergence error
	return (u[p] - uk) * (u[p] - uk) + (v[p] - vk) * (v[p] - vk);
}



/**
 *
 *  Horn & Schunck method for optical flow estimation at a single scale
 *
 */
void horn_schunck_optical_flow(
	const float *I1,             // source image
	const float *I2,             // target image
	float       *u,              // x component of optical flow
	float       *v,              // y component of optical flow
	const int    nx,             // image width
	const int    ny,             // image height
	const float  alpha,          // smoothing parameter
	const int    warps,          // number of warpings per scale
	const float  TOL,            // stopping criterion threshold
	const int    maxiter,        // maximum number of iterations
	const bool   verbose         // switch on messages
)
{
	if (verbose) fprintf(stderr, "Single-scale Horn-Schunck of a %dx%d "
		       "image\n\ta=%g nw=%d eps=%g mi=%d v=%d\n", nx, ny,
			alpha, warps, TOL, maxiter, verbose);

	const int   size   = nx * ny;
	const float alpha2 = alpha * alpha;

	//allocate memory
	int sf = sizeof(float);
	float *I2x  = xmalloc(size * sf); // x derivative of I2
	float *I2y  = xmalloc(size * sf); // y derivative of I2
	float *I2w  = xmalloc(size * sf); // warping of I2
	float *I2wx = xmalloc(size * sf); // warping of I2x
	float *I2wy = xmalloc(size * sf); // warping of I2y
	float *Au   = xmalloc(size * sf); // constant part of numerator of u
	float *Av   = xmalloc(size * sf); // constant part of numerator of v
	float *Du   = xmalloc(size * sf); // denominator of u
	float *Dv   = xmalloc(size * sf); // denominator of v
	float *D    = xmalloc(size * sf); // common numerator of u and v

	// compute the gradient of the second image
	gradient(I2, I2x, I2y, nx, ny);

	// iterative approximation to the Taylor expansions
	for(int n = 0; n < warps; n++)
	{
		if(verbose) fprintf(stderr, "Warping %d:", n);

		// warp the second image and its derivatives
		bicubic_interpolation_warp(I2,  u, v, I2w,  nx, ny, true);
		bicubic_interpolation_warp(I2x, u, v, I2wx, nx, ny, true);
		bicubic_interpolation_warp(I2y, u, v, I2wy, nx, ny, true);

		// store the constant parts of the system
		for(int i = 0; i < size; i++)
		{
			const float I2wl = I2wx[i] * u[i] + I2wy[i] * v[i];
			const float dif  = I1[i] - I2w[i] + I2wl;

			Au[i] = dif * I2wx[i];
			Av[i] = dif * I2wy[i];
			Du[i] = I2wx[i] * I2wx[i] + alpha2;
			Dv[i] = I2wy[i] * I2wy[i] + alpha2;
			D[i]  = I2wx[i] * I2wy[i];
		}

		int niter = 0;
		float error = 1000;

		// iterations of the SOR numerical scheme
		while(error > TOL && niter < maxiter)
		{
			niter++;
			error = 0;

			//process the central part of the optical flow
			#pragma omp parallel for reduction(+:error)
			for(int i = 1; i < ny-1; i++)
			for(int j = 1; j < nx-1; j++)
			{
				const int k = i * nx + j;
				error += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx-1, k-nx+1, k+nx-1,
						k+nx+1, k-nx, k-1, k+nx, k+1
						);
			}

			// process the first and last rows
			for(int j = 1; j < nx-1; j++)
			{
				// first row
				int k = j;
				error += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-1, k+1, k+nx-1, k+nx+1,
						k, k-1, k+nx, k+1
						);

				// last row
				k = (ny-1) * nx + j;
				error += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx-1, k-nx+1, k-1, k+1,
						k-nx, k-1, k, k+1
						);
			}

			// process the first and last columns
			for(int i = 1; i < ny-1; i++)
			{
				// first column
				int k = i * nx;
				error += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx, k-nx+1, k+nx, k+nx+1,
						k-nx, k, k+nx, k+1
						);

				// last column
				k = (i+1) * nx - 1;
				error += sor_iteration(
						Au, Av, Du, Dv, D, u, v, alpha2,
						k, k-nx-1, k-nx, k+nx-1, k+nx,
						k-nx, k-1, k+nx, k
						);
			}

			// process the corners
			// up-left corner
			error += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					0, 0, 1, nx, nx+1,
					0, 0, nx, 1
					);

			// up-right corner
			int k = nx - 1;
			error += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					k, k-1, k, k+nx-1, k+nx,
					k, k-1, k+nx, k
					);

			// bottom-left corner
			k = (ny-1) * nx;
			error += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					k, k-nx, k-nx+1,k, k+1,
					k-nx, k, k, k+1
					);

			// bottom-right corner
			k = ny * nx - 1;
			error += sor_iteration(
					Au, Av, Du, Dv, D, u, v, alpha2,
					k, k-1, k, k-nx-1, k-nx,
					k-nx, k-1, k, k
					);

			error = sqrt(error / size);
		}

		if(verbose)
			fprintf(stderr, "Iterations %d (%g)\n", niter, error);
	}

	// free the allocated memory
	free(I2x);
	free(I2y);
	free(I2w);
	free(I2wx);
	free(I2wy);
	free(Au);
	free(Av);
	free(Du);
	free(Dv);
	free(D);
}

// compute the largest number of an array
static float max_element(const float *x, int n)
{
	int r = 0;
	for (int i = 1; i < n; i++)
		if (x[i] > x[r])
			r = i;
	return x[r];
}

// compute the smallest number of an array
static float min_element(const float *x, int n)
{
	int r = 0;
	for (int i = 1; i < n; i++)
		if (x[i] < x[r])
			r = i;
	return x[r];
}


/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
static void image_normalization(
	const float *I1,   // first image
	const float *I2,   // second image
	float       *I1n,  // first normalized image
	float       *I2n,  // second normalized image
	int          size  // size of each image
)
{
	// find the max and min of both images
	const float max1 = max_element(I1, size);
	const float max2 = max_element(I2, size);
	const float min1 = min_element(I1, size);
	const float min2 = min_element(I2, size);

	// obtain the absolute max and min
	const float max = max1 > max2 ? max1 : max2;
	const float min = min1 < min2 ? min1 : min2;
	const float den = max - min;

	if(den > 0)
		// normalize both images
		for(int i = 0; i < size; i++)
		{
			I1n[i] = 255.0 * (I1[i] - min) / den;
			I2n[i] = 255.0 * (I2[i] - min) / den;
		}

	else
		// copy the original images
		for(int i = 0; i < size; i++)
		{
			I1n[i] = I1[i];
			I2n[i] = I2[i];
		}
}


/**
 *
 *  Procedure to handle the pyramidal approach.
 *  This procedure relies on the previous functions to calculate
 *  large optical flow fields using a pyramidal scheme.
 *
 */
 void horn_schunck_pyramidal(
	const float *I1,              // source image
	const float *I2,              // target image
	float       *u,               // x component of optical flow
	float       *v,               // y component of optical flow
	const int    nx,              // image width
	const int    ny,              // image height
	const float  alpha,           // smoothing weight
	const int    nscales,         // number of scales
	const float  zfactor,         // zoom factor
	const int    warps,           // number of warpings per scale
	const float  TOL,             // stopping criterion threshold
	const int    maxiter,         // maximum number of iterations
	const bool   verbose          // switch on messages
)
{
	if (verbose) fprintf(stderr, "Multiscale Horn-Schunck of a %dx%d pair"
			"\n\ta=%g ns=%d zf=%g nw=%d eps=%g mi=%d\n", nx, ny,
			alpha, nscales, zfactor, warps, TOL, maxiter);

	int size = nx * ny;

	float *I1s[nscales];
	float *I2s[nscales];
	float *us[nscales];
	float *vs[nscales];
	int nxx[nscales];
	int nyy[nscales];


	I1s[0] = xmalloc(size * sizeof(float));
	I2s[0] = xmalloc(size * sizeof(float));

	// normalize the finest scale images between 0 and 255
	image_normalization(I1, I2, I1s[0], I2s[0], size);

	// presmoothing the finest scale images
	gaussian(I1s[0], nx, ny, INPUT_PRESMOOTHING_SIGMA);
	gaussian(I2s[0], nx, ny, INPUT_PRESMOOTHING_SIGMA);

	us[0] = u;
	vs[0] = v;
	nxx[0] = nx;
	nyy[0] = ny;

	// create the scales
	for(int s = 1; s < nscales; s++)
	{
		zoom_size(nxx[s-1], nyy[s-1], nxx+s, nyy+s, zfactor);
		const int sizes = nxx[s] * nyy[s];

		I1s[s] = xmalloc(sizes * sizeof(float));
		I2s[s] = xmalloc(sizes * sizeof(float));
		us[s] = xmalloc(sizes * sizeof(float));
		vs[s] = xmalloc(sizes * sizeof(float));

		// compute the zoom from the previous finer scale
		zoom_out(I1s[s-1], I1s[s], nxx[s-1], nyy[s-1], zfactor);
		zoom_out(I2s[s-1], I2s[s], nxx[s-1], nyy[s-1], zfactor);
	}

	// initialize the flow
	for (int i = 0; i < nxx[nscales-1] * nyy[nscales-1]; i++)
	{
		us[nscales-1][i] = 0;
		vs[nscales-1][i] = 0;
	}

	// pyramidal approximation to the optic flow
	for(int s = nscales-1; s >= 0; s--)
	{
		if(verbose)
			fprintf(stderr, "Scale: %d %dx%d\n", s, nxx[s], nyy[s]);

		// compute the optical flow at this scale
		horn_schunck_optical_flow(
			I1s[s], I2s[s], us[s], vs[s], nxx[s], nyy[s],
			alpha, warps, TOL, maxiter, verbose
		);

		// if this was the last scale, finish now
		if (!s) break;

		// otherwise, upsample the optical flow

		// zoom the optic flow for the next finer scale
		zoom_in(us[s], us[s-1], nxx[s], nyy[s], nxx[s-1], nyy[s-1]);
		zoom_in(vs[s], vs[s-1], nxx[s], nyy[s], nxx[s-1], nyy[s-1]);

		// scale the optic flow with the appropriate zoom factor
		for(int i = 0; i < nxx[s-1] * nyy[s-1]; i++)
		{
			us[s-1][i] *= 1.0 / zfactor;
			vs[s-1][i] *= 1.0 / zfactor;
		}
	}

	// free the allocated memory
	free(I1s[0]);
	free(I2s[0]);
	for(int i = 1; i < nscales; i++)
	{
		free(I1s[i]);
		free(I2s[i]);
		free(us[i]);
		free(vs[i]);
	}
}


#endif//HORN_SCHUNCK_PYRAMIDAL_C
