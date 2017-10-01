// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
// <jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
// All rights reserved.

#ifndef CLG_OF_H
#define CLG_OF_H

#define EPS 1E-12           // Precision threshold for linear system solving.
#define MIN_ERROR 0.0001    // Threshold to stop iterations, checked against
                            // the mean variation of the flow field.
#define MIN_GAUSSIAN_SIZE 3 // Minimum filter size required to perform
                            // gaussian smoothing of the input images.

#include "util.h"
#include "zoom.h"
#include "iio.h"

// Neumann boundary condition (derivatives are set to zero).
double boundaryCondition(double **u, double **v, int nRows, int nCols);

double SOR_at(double **u,
              double **v,
              double **J[JROWS][JCOLS],
              int i,
              int j,
              double alpha,
              double wFactor);

// Coupled Gauss-Seidel iteration u^{k+1} and v^{k+1} are simultaneously
// computed, improving convergence speed.
double relaxPointwiseCoupledGaussSeidel(double **u,
                                        double **v,
                                        double **J[JROWS][JCOLS],
                                        int nRows,
                                        int nCols,
                                        double alpha);

// SOR iteration whre u^{k+1} and v^{k+1} are sequentially computed.
// If w=1.0 this method is equivalent to standard GaussSeidel.
double relaxSOR(double **u,
                double **v,
                double **J[JROWS][JCOLS],
                int nRows,
                int nCols,
                double alpha,
                double wFactor);

// Optical flow computation, with numIteration calls to either SOR
// or Pointwise-Coupled Gauss-Seidel.
int calcCLG_OF(double* image1,
               double* image2,
               double* uOut,
               double* vOut,
               int nRows,
               int nCols,
               int iterations,
               double alpha,
               double rho,
               double wFactor,
               int verbose,
               int coupledMode);

// Multiscale optical flow. For each nScale, there is a call of the
// calcCLG_OF function. Warping is performed by bicubic interpolation.
int calcMSCLG_OF(double* image1,
                 double* image2,
                 double* uOut,
                 double* vOut,
                 int nRows,
                 int nCols,
                 int iterations,
                 double alpha,
                 double rho,
                 double sigma,
                 double wFactor,
                 int nScales,
                 const double scaleFactor,
                 int coupledMode,
                 int verbose);

#endif
