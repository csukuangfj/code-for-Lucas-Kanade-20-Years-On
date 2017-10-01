// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef ZOOM_H
#define ZOOM_H

#define ZOOM_SIGMA_ZERO 0.6

#include "util.h"
#include "mask.h"
#include "bicubic_interpolation.h"

void zoom_size(
    int nx,             // width of the orignal image
    int ny,             // height of the orignal image
    int *nxx,           // width of the zoomed image
    int *nyy,           // height of the zoomed image
    double factor        // zoom factor between 0 and 1
);

void zoom_out(
    const double *I,          // input image
    double *Iout,             // output image
    const int nx,            // image width
    const int ny,            // image height
    const double factor       // zoom factor between 0 and 1
);

void zoom_in(
    const double *I, // input image
    double *Iout,    // output image
    int nx,         // width of the original image
    int ny,         // height of the original image
    int nxx,        // width of the zoomed image
    int nyy         // height of the zoomed image
);



#endif//ZOOM_H
