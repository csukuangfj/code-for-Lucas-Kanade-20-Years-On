// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright 2012, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>
// All rights reserved.



#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// this function is like "malloc", but it returns always a valid pointer
static void *xmalloc(size_t size)
{
	void *new = malloc(size);
	if (!new)
		exit(fprintf(stderr, "out of memory\n"));
	return new;
}

// the type of the "getpixel" function
typedef float (*extension_operator_float)(float*, int, int, int, int);

// getpixel, with neumann boundary conditions
static float extend_float_image_constant(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

// compute the gradient and temporal derivative of the input image pair
static void compute_input_derivatives(float *Ex, float *Ey, float *Et,
		float *a, float *b, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++) {
		Ey[j*w+i] = (1.0/4) * ( p(a,w,h, i, j+1) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i+1,j));
		Ex[j*w+i] = (1.0/4) * ( p(a,w,h, i+1, j) - p(a,w,h, i,j)
				+ p(a,w,h, i+1, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j) - p(b,w,h, i,j)
				+ p(b,w,h, i+1, j+1) - p(b,w,h, i,j+1));
		Et[j*w+i] = (1.0/4) * ( p(b,w,h, i, j) - p(a,w,h, i,j)
				+ p(b,w,h, i+1, j) - p(a,w,h, i+1,j)
				+ p(b,w,h, i, j+1) - p(a,w,h, i,j+1)
				+ p(b,w,h, i+1, j+1) - p(a,w,h, i+1,j+1));
	}
}

// compute a local average of a function "u"
static void compute_bar(float *ubar, float *u, int w, int h)
{
	extension_operator_float p = extend_float_image_constant;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		ubar[j*w+i] = (1.0/6) * (p(u,w,h, i-1, j) + p(u,w,h, i+1, j)
				+ p(u,w,h, i, j-1) + p(u,w,h, i, j+1))
			+ (1.0/12) * (p(u,w,h, i-1,j-1) + p(u,w,h, i+1,j-1)
				+ p(u,w,h, i-1,j+1) + p(u,w,h, i+1,j+1));
}

// compute a sigle iteration of the classical Horn-Schunck method
static void hs_iteration(float *u, float *v,
		float *Ex, float *Ey, float *Et, int w, int h, float alpha)
{
	float *ubar = xmalloc(w * h * sizeof(float));
	float *vbar = xmalloc(w * h * sizeof(float));
	compute_bar(ubar, u, w, h);
	compute_bar(vbar, v, w, h);
	for (int i = 0; i < w*h; i++) {
		float t = Ex[i]*ubar[i] + Ey[i]*vbar[i] + Et[i];
		t /= alpha*alpha + Ex[i]*Ex[i] + Ey[i]*Ey[i];
		u[i] = ubar[i] - Ex[i] * t;
		v[i] = vbar[i] - Ey[i] * t;
	}
	free(ubar);
	free(vbar);
}

// run n iterations of the classical Horn-Schunck method
void hs(float *u, float *v, float *a, float *b, int w, int h,
		int n, float alpha)
{
	float *gx = xmalloc(w * h * sizeof(float));
	float *gy = xmalloc(w * h * sizeof(float));
	float *gt = xmalloc(w * h * sizeof(float));
	compute_input_derivatives(gx, gy, gt, a, b, w, h);
	for (int i = 0; i < w*h; i++)
		u[i] = v[i] = 0;
	for (int i = 0; i < n; i++)
		hs_iteration(u, v, gx, gy, gt, w, h, alpha);
	free(gx);
	free(gy);
       	free(gt);
}


#ifndef OMIT_MAIN
#include "iio.h"

// main function for testing the Horn-Schunck method from the command line
int main(int argc, char *argv[])
{
	if (argc != 6 && argc != 7)
		return fprintf(stderr,"usage:\n\t%s niter alpha a b f\n",*argv);
	//                                        0 1     2     3 4 5
	int niter = atoi(argv[1]);
	float alpha = atof(argv[2]);
	char *filename_a = argv[3];
	char *filename_b = argv[4];
	char *filename_f = argv[5];
	int w, h, ww, hh;
	float *a = iio_read_image_float(filename_a, &w, &h);
	float *b = iio_read_image_float(filename_b, &ww, &hh);
	if (w != ww || h != hh)
		return fprintf(stderr, "input images size mismatch\n");
	float *u = xmalloc(2 * w * h * sizeof(float));
	float *v = u + w*h;
	hs(u, v, a, b, w, h, niter, alpha);
	iio_save_image_float_split(filename_f, u, w, h, 2);
	return EXIT_SUCCESS;
}
#endif//OMIT_MAIN
