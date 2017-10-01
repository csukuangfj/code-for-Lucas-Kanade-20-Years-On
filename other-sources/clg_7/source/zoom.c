// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef ZOOM_C
#define ZOOM_C

#include "zoom.h"

/**
  *
  * Compute the size of a zoomed image from the zoom factor
  *
**/
void zoom_size(
    int nx,             // width of the orignal image
    int ny,             // height of the orignal image
    int *nxx,           // width of the zoomed image
    int *nyy,           // height of the zoomed image
    double factor        // zoom factor between 0 and 1
)
{
    *nxx = (int)((double) nx * factor + 0.5);
    *nyy = (int)((double) ny * factor + 0.5);
}

/**
  *
  * Downsample an image
  *
**/
void zoom_out(
    const double *I,          // input image
    double *Iout,             // output image
    const int nx,            // image width
    const int ny,            // image height
    const double factor       // zoom factor between 0 and 1
)
{
    // temporary working image
    double *Is = xmalloc(nx * ny * sizeof*Is);
    for(int i = 0; i < nx * ny; i++)
        Is[i] = I[i];

    // compute the size of the zoomed image
    int nxx, nyy;
    zoom_size(nx, ny, &nxx, &nyy, factor);

    // compute the Gaussian sigma for smoothing
    const double sigma = ZOOM_SIGMA_ZERO * sqrt(1.0/(factor*factor) - 1.0);

    // pre-smooth the image
    gaussian(Is, nx, ny, sigma);

    // re-sample the image using bicubic interpolation
    //#pragma omp parallel for
    for (int i1 = 0; i1 < nyy; i1++)
    for (int jj1 = 0; jj1 < nxx; jj1++)
    {
        const double i2  = (double) i1 / factor;
        const double j2  = (double) jj1 / factor;

        double g = bicubic_interpolation_at(Is, j2, i2, nx, ny, false);
        Iout[i1 * nxx + jj1] = g;
    }

    free(Is);
}


/**
  *
  * Function to upsample the image
  *
**/
void zoom_in(
    const double *I, // input image
    double *Iout,    // output image
    int nx,         // width of the original image
    int ny,         // height of the original image
    int nxx,        // width of the zoomed image
    int nyy         // height of the zoomed image
)
{
    // compute the zoom factor
    const double factorx = ((double)nxx / nx);
    const double factory = ((double)nyy / ny);

    // re-sample the image using bicubic interpolation
    //#pragma omp parallel for
    for (int i1 = 0; i1 < nyy; i1++)
    for (int jj1 = 0; jj1 < nxx; jj1++)
    {
        double i2 =  (double) i1 / factory;
        double j2 =  (double) jj1 / factorx;

        double g = bicubic_interpolation_at(I, j2, i2, nx, ny, false);
        Iout[i1 * nxx + jj1] = g;
    }
}



#endif//ZOOM_C
