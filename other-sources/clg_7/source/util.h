// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
// <jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
// All rights reserved.

#ifndef UTIL_H
#define UTIL_H

#define JROWS 3     // Array size for storing J coefficient values.
#define JCOLS 3

#include <stdio.h>  // Used for "printf".
#include <stdlib.h> // Used for "free".
#include <math.h>   // Used for "exp".
#include "mask.h"   // Used for 2D gaussian convolution.
#include "zoom.h"

// 2x2 matrix determinant.
double lin2by2det(double a, double b, double c, double d);

// Square matrix memory allocation.
double **pMatrix(int nRows, int nCols);

// Square matrix memory free.
void freePmatrix(double **mat, int nRows);

// Verifies % correct if an index is valid in the range [0,size-1].
int correctIndex(int p, int size);

// Derivatives for each pixel.
void computeDerivatives(double **image1,
                        double **image2,
                        double **dfdx,
                        double **dfdy,
                        int nRows,
                        int nCols,
                        int verbose);

// Sets the J tensor values.
void computeJTensor(double **dfdx,
                    double **dfdy,
                    double **dfdt,
                    double **J[JROWS][JCOLS],
                    int nRows,
                    int nCols);

// In-place gaussian convolution by using mask.h implementation.
void matrixSmooth(double **matrix,
                  int nRows,
                  int nCols,
                  double kernelSigma);

// Array memory allocation.
void *xmalloc(size_t size);

#endif
