// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef MASK_H
#define MASK_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

#define BOUNDARY_CONDITION_DIRICHLET 0
#define BOUNDARY_CONDITION_REFLECTING 1
#define BOUNDARY_CONDITION_PERIODIC 2

#define DEFAULT_GAUSSIAN_WINDOW_SIZE 5.0
#define DEFAULT_BOUNDARY_CONDITION BOUNDARY_CONDITION_REFLECTING


void gradient(
    const double *input, // input image
    double *dx,          // computed x derivative
    double *dy,          // computed y derivative
    const int nx,       // image width
    const int ny);        // image height

void gaussian(
    double *I,            // input/output image
    const int xdim,       // image width
    const int ydim,       // image height
    const double sigma    // Gaussian sigma
);


#endif//MASK_C
