// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
// <jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
// All rights reserved.

#ifndef CLG_OF_C
#define CLG_OF_C

/**
 * @file clg_of.c
 * @brief Implementation of CLG optical flow for 2D grayscale images.
 * @author Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano <jdelpian@uandes.cl>
 */
#include "clg_of.h"
#include <math.h>   // Used for "exp".
#include <stdio.h>  // Used for "printf".
#include <stdlib.h> // Used for "free".


/**
 * boundaryCondition
 *
 * Neumann boundary conditions (derivatives are set to zero).
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 *
 */
double boundaryCondition(double **u, double **v, int nRows, int nCols) {

    int i, j;
    double error = 0.0;

    // First and last rows.
    for (j=1; j<nCols-1; j++) {

        i = 0;
        error += (u[i+1][j]-u[i][j]) * (u[i+1][j]-u[i][j]);
        error += (v[i+1][j]-v[i][j]) * (v[i+1][j]-v[i][j]);

        u[i][j] = u[i+1][j];
        v[i][j] = v[i+1][j];

        i = nRows-1;
        error += (u[i-1][j]-u[i][j]) * (u[i-1][j]-u[i][j]);
        error += (v[i-1][j]-v[i][j]) * (v[i-1][j]-v[i][j]);

        u[i][j] = u[i-1][j];
        v[i][j] = v[i-1][j];
    }

    // First and last columns.
    for (i=1; i<nRows-1; i++) {

        j = 0;
        error += (u[i][j+1]-u[i][j]) * (u[i][j+1]-u[i][j]);
        error += (v[i][j+1]-v[i][j]) * (v[i][j+1]-v[i][j]);

        u[i][j] = u[i][j+1];
        v[i][j] = v[i][j+1];

        j = nCols-1;
        error += (u[i][j-1]-u[i][j]) * (u[i][j-1]-u[i][j]);
        error += (v[i][j-1]-v[i][j]) * (v[i][j-1]-v[i][j]);

        u[i][j] = u[i][j-1];
        v[i][j] = v[i][j-1];
    }

    // Corners.
    i=0; j=0;
    error += (u[i+1][j+1]-u[i][j]) * (u[i+1][j+1]-u[i][j]);
    error += (v[i+1][j+1]-v[i][j]) * (v[i+1][j+1]-v[i][j]);

    u[i][j] = u[i+1][j+1];
    v[i][j] = v[i+1][j+1];

    j = nCols-1;
    error += (u[i+1][j-1]-u[i][j]) * (u[i+1][j-1]-u[i][j]);
    error += (v[i+1][j-1]-v[i][j]) * (v[i+1][j-1]-v[i][j]);

    u[i][j] = u[i+1][j-1];
    v[i][j] = v[i+1][j-1];

    i=nRows-1; j=0;
    error += (u[i-1][j+1]-u[i][j]) * (u[i-1][j+1]-u[i][j]);
    error += (v[i-1][j+1]-v[i][j]) * (v[i-1][j+1]-v[i][j]);

    u[i][j] = u[i-1][j+1];
    v[i][j] = v[i-1][j+1];

    j = nCols-1;
    error += (u[i-1][j-1]-u[i][j]) * (u[i-1][j-1]-u[i][j]);
    error += (v[i-1][j-1]-v[i][j]) * (v[i-1][j-1]-v[i][j]);

    u[i][j] = u[i-1][j-1];
    v[i][j] = v[i-1][j-1];

    return error;
}


/**
 * SOR_at
 *
 * SOR iteration at location (i,j).
 * This method return the new values.
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
 *                 smoothed) derivatives for each image/optical flow pixel.
 * i               Row of the location to compute.
 * j               Column of the location to compute.
 * alpha           Optical flow global smoothing coefficient.
 * wFactor         Relaxation parameter (if w=1.0 then the function becomes the
 *                 Gauss-Seidel method).
 *
 */
double SOR_at(double **u,
              double **v,
              double **J[JROWS][JCOLS],
              int i,
              int j,
              double alpha,
              double wFactor) {

    double h2a, numerator, denominator, ulocal, vlocal, error;

    // h, could be less than 1.0.
    h2a = 1.0 / alpha;

    // SOR formula.
    numerator   = u[i][j-1] + u[i-1][j] + u[i][j+1] + u[i+1][j] -
                  h2a * (J[0][1][i][j] * v[i][j]+J[0][2][i][j]);

    denominator = 4.0 + h2a * J[0][0][i][j];

    ulocal      = (1.0-wFactor) * u[i][j] + wFactor * numerator / denominator;

    numerator   = v[i][j-1] + v[i-1][j] + v[i][j+1] + v[i+1][j] -
                  h2a * (J[0][1][i][j] * ulocal+J[1][2][i][j]);

    denominator = 4.0 + h2a * J[1][1][i][j];

    vlocal      = (1.0-wFactor) * v[i][j] + wFactor * numerator / denominator;

    error       = (ulocal - u[i][j]) * (ulocal - u[i][j]);
    error      += (vlocal - v[i][j]) * (vlocal - v[i][j]);

    u[i][j]     = ulocal;
    v[i][j]     = vlocal;

    return error;
}


/**
 * relaxPointwiseCoupledGaussSeidel
 *
 * Pointwise coupled Gauss-Seidel relaxation iteration for CLG-OF equations.
 * Each call to this function updates the current value of the solution,
 * u[1..m][1..n], v[1..m][1..n], using the motion tensor
 * J[JROWS][JCOLS][1..m][1..n].
 * Neumann boundary conditions are used (derivatives are set to zero).
 * The return value is the total error of the current iteration.
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
 *                 smoothed) derivatives for each image/optical flow pixel.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 * alpha           Optical flow global smoothing coefficient.
 *
 */
double relaxPointwiseCoupledGaussSeidel(double **u,
                                        double **v,
                                        double **J[JROWS][JCOLS],
                                        int nRows,
                                        int nCols,
                                        double alpha) {
    int i, j;
    double hx, hy, hx2, hy2, discr, detU, detV, error, ulocal, vlocal;

    // h, could be less than 1.0.
    hx    = 1.0;
    hy    = hx;
    hx2   = hx * hx;
    hy2   = hy * hy;
    error = 0.0;

    for (i=1; i<nRows-1; i++) {
        for (j=1; j<nCols-1; j++) {
            // Point-wise Coupled Gauss-Seidel formula.

            // Coupled system for two variables solved by Cramer's rule.
            // This is the discriminant of that system
            discr = lin2by2det(alpha/hx2*2.0 + alpha/hy2*2.0 + J[0][0][i][j],
                               J[0][1][i][j],
                               J[0][1][i][j],
                               alpha/hx2*2.0 + alpha/hy2*2.0 + J[1][1][i][j]);

            if (abs(discr) > EPS) {
                // Determinant replacing first column to solve
                // the 2x2 system by using Cramer's rule.
                detU = lin2by2det(alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                  + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                  - J[0][2][i][j],
                                    J[0][1][i][j],
                                    alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                  + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                  - J[1][2][i][j],
                                    alpha / hx2 * 2.0 + alpha / hy2 * 2.0
                                  + J[1][1][i][j]);

                // Determinant replacing second column to solve
                // the 2x2 system by using Cramer's rule.
                detV = lin2by2det(alpha / hx2 * 2.0 + alpha / hy2 * 2.0
                                  + J[0][0][i][j],
                                    alpha / hx2 * (u[i-1][j] + u[i+1][j])
                                  + alpha / hy2 * (u[i][j-1] + u[i][j+1])
                                  - J[0][2][i][j],
                                    J[0][1][i][j],
                                    alpha / hx2 * (v[i-1][j] + v[i+1][j])
                                  + alpha / hy2 * (v[i][j-1] + v[i][j+1])
                                  - J[1][2][i][j]);

                // Division of two discriminants (Cramer's rule).
                ulocal = detU / discr;
                vlocal = detV / discr;

                error += (ulocal-u[i][j])*(ulocal-u[i][j]);
                error += (vlocal-v[i][j])*(vlocal-v[i][j]);

                u[i][j]=ulocal;
                v[i][j]=vlocal;

                // Normal Gauss-Seidel iteration.
            } else {

                // w_Factor=1.0, thus a Gauss-Seidel iteration.
                error += SOR_at(u, v, J, i, j, alpha, 1.0);

            } // End if.
        } // End columns.
    } // End rows.

    error += boundaryCondition(u, v, nRows, nCols);

    return sqrt(error / (nRows*nCols));;
}


/**
 * relaxSOR
 *
 * SOR relaxation iteration for CLG-OF equations.
 * Each call to this function updates the current value of the solution,
 * u[1..m][1..n], v[1..m][1..n], using the motion tensor
 * J[JROWS][JCOLS][1..m][1..n].
 * Neumann boundary conditions are used (derivatives are set to zero).
 *
 * Parameters:
 *
 * u               Pointer to the optical flow horizontal component
 *                 matrix/array.
 * v               Pointer to the optical flow vertical component matrix/array.
 * J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
 *                 smoothed) derivatives for each image/optical flow pixel.
 * nRows           Number of rows of the optical flow arrrays.
 * nCols           Number of columns of the optical flow arays.
 * alpha           Optical flow global smoothing coefficient.
 * wFactor         Relaxation parameter (if w=1.0 then the function is the
 *                 Gauss-Seidel method).
 *
 */
double relaxSOR(double **u,
                double **v,
                double **J[JROWS][JCOLS],
                int nRows,
                int nCols,
                double alpha,
                double wFactor) {
    int i, j;
    double error = 0.0;

    for (i=1; i<nRows-1; i++)
        for (j=1; j<nCols-1; j++)
            error += SOR_at(u, v, J, i, j, alpha, wFactor);

    error += boundaryCondition(u, v, nRows, nCols);

    return sqrt(error / (nRows*nCols));
}


/**
 * calcCLG_OF
 *
 * Main CLG-optical flow (CLG-OF) computation function.
 *
 * Parameters:
 *
 * image1          Pointer to the first image (the "previous" time frame).
 * image2          Pointer to the second image (the "current" time frame).
 * uOut            Pointer to the horizontal component of the CLG-OF solution.
 * vOut            Pointer to the vertical component of the CLG-OF solution.
 * nRows           Number of image rows (same for the CLG-OF vector field).
 * nCols           Number of image columns (same for the CLG-OF vector field).
 * iterations      Number of iterations for iterative solution.
 * alpha           Global smoothing coefficient of the CLG-OF.
 * rho             Local spatio-temporal smoothing coefficient of the CLG-OF.
 * wFactor         SOR relaxation factor, between 0 and 2.
 * verbose         Display/hide messages to the stdout.
 * coupledMode     Iteration type. 1->Pointwise-Coupled Gauss-Seidel, 0->SOR.
 *
 */
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
               int coupledMode) {

    if (verbose) {
        printf("calc_clg\n");
        printf("  setting up variables\n");
    }

    int i=0, j=0;

    // h, could be less than 1.0.
    double h = 1.0;

    // Matrix to vector.
    double **prevFrame, **currFrame;
    prevFrame = pMatrix(nRows, nCols);
    currFrame = pMatrix(nRows, nCols);

    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            prevFrame[i][j] = image1[j + i*nCols];
            currFrame[i][j] = image2[j + i*nCols];
        }
    }

    if (verbose)
        printf("  allocating memory for arrays\n");

    double **u, **v;
    double **dfdx, **dfdxw;
    double **dfdy, **dfdyw;
    double **dfdt, **dfdtw;

    u = pMatrix(nRows, nCols);
    v = pMatrix(nRows, nCols);

    // Derivatives and their warped versions.
    dfdx = pMatrix(nRows, nCols);
    dfdy = pMatrix(nRows, nCols);
    dfdt = pMatrix(nRows, nCols);

    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            dfdt[i][j] = currFrame[i][j] - prevFrame[i][j];
        }
    }

    if (verbose)
        printf("  allocating memory for derivatives matrices\n");

    double **J[JROWS][JCOLS];

    // Because of symmetry, only the upper part is allocated.
    int k, l;
    for (k=0; k<JROWS; k++)
        for (l=k; l<JCOLS; l++)
            J[k][l] = pMatrix(nRows, nCols);

    // Spatial derivatives obtention.
    computeDerivatives(prevFrame, currFrame, dfdx, dfdy, nRows, nCols, verbose);
    // Compute J tensor.
    computeJTensor(dfdx, dfdy, dfdt, J, nRows, nCols);

    if (verbose)
        printf("  local spatio temporal smoothing\n");

    if (rho > 0) {

        k=0, l=0;
        for (k=0; k<JROWS; k++)
            for (l=k; l<JCOLS; l++)
                matrixSmooth(J[k][l], nRows, nCols, rho);

    }

    if (iterations == 0)
        iterations = (int) (nRows * nCols / 8.0);

    if (verbose)
        printf("  performing %i relax iterations\n", iterations);

    int count = 0;
    double error = 1000000;
    double convergenceError = 0.0;

    for (count=0; count<iterations && convergenceError*1.01 < error; count++) {

        if (count > 0)
            error = convergenceError;

        if (coupledMode == 1) {

            if (verbose && count % 50 == 0 && count > 0)
                printf("  iteration %d/%d (P-C Gauss-Seidel), error=%f\n",
                       count, iterations, error);

            convergenceError = relaxPointwiseCoupledGaussSeidel(u, v, J,
                                                                nRows, nCols,
                                                                alpha);
        } else {
            if (verbose && count % 50 == 0 && count > 0)
                printf("  iteration %d/%d (SOR), error=%f\n",
                       count, iterations, error);

            convergenceError = relaxSOR(u, v, J, nRows, nCols, alpha, wFactor);
        }
    }

    // Show debug information.
    if (verbose)
        printf("  filling output after %d iterations, error=%f\n", count, error);

    // Fill output variables.
    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            uOut[j + i*nCols] += u[i][j];
            vOut[j + i*nCols] += v[i][j];
        }
    }

    // Free memory.
    if (verbose)
        printf("  freeing memory\n");

    freePmatrix(u, nRows);
    freePmatrix(v, nRows);

    for (k=0; k<JROWS; k++)
        for (l=k; l<JCOLS; l++)
            freePmatrix(J[k][l], nRows);

    freePmatrix(prevFrame, nRows);
    freePmatrix(currFrame, nRows);
    freePmatrix(dfdx,  nRows);
    freePmatrix(dfdy,  nRows);
    freePmatrix(dfdt,  nRows);

    if (verbose)
        printf("calc_clg: done\n");
    return 1;
}


/**
 * calcMSCLG_OF
 *
 * Main Multiscale CLG-optical flow (CLG-OF) computation function.
 *
 * Parameters:
 *
 * image1          Pointer to the first image (the "previous" time frame).
 * image2          Pointer to the second image (the "current" time frame).
 * uOut            Pointer to the horizontal component of the CLG-OF solution.
 * vOut            Pointer to the vertical component of the CLG-OF solution.
 * nRows           Number of image rows (same for the CLG-OF vector field).
 * nCols           Number of image columns (same for the CLG-OF vector field).
 * iterations      Number of iterations for iterative solution.
 * alpha           Global smoothing coefficient of the CLG-OF.
 * rho             Local spatio-temporal smoothing coefficient of the CLG-OF.
 * sigma           Standard deviation for an optional gaussian smoothing kernel,
 *                 applied to the input images prior to the CLG-OF computation.
 * wFactor         SOR relaxation factor, between 0 and 2.
 * nScales         Total number of scales (at least 1).
 * scaleFactor     Downsampling factor, between 0.5 and 1.0.
 * coupledMode     Iteration type. 1->Pointwise-Coupled Gauss-Seidel, 0->SOR.
 * verbose         Display/hide messages to the stdout.
 *
 */
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
                 int verbose) {

    if (verbose) {
        printf("calcMSCLG_OF\n");
        printf("  parameters: nScales=%d, scaleFactor=%f\n",
               nScales, scaleFactor);
    }

    // Variables declaration.
    int s, i;
    double *image2warped;
    int size = nCols * nRows;

    double *I1s[nScales];
    double *I2s[nScales];
    double *us[nScales];
    double *vs[nScales];
    int nxx[nScales];
    int nyy[nScales];

    I1s[0] = image1;
    I2s[0] = image2;

    // Pre-smoothing the finest scale images.
    int filterSize = (int) 2*(DEFAULT_GAUSSIAN_WINDOW_SIZE * sigma) + 1;

    if (filterSize >= MIN_GAUSSIAN_SIZE) {
        printf("  apply gaussian smoothing for each frame\n");
        gaussian(I1s[0], nCols, nRows, sigma);
        gaussian(I2s[0], nCols, nRows, sigma);
    } else {
        printf("  no gaussian smoothing was applied to input frames\n");
    }

    us[0] = uOut;
    vs[0] = vOut;
    nxx[0] = nCols;
    nyy[0] = nRows;

    // Create the scales.
    for (s=1; s<nScales; s++) {

        if (verbose)
            printf("  scales %d init\n", s);

        zoom_size(nxx[s-1], nyy[s-1], nxx+s, nyy+s, scaleFactor);
        const int sizes = nxx[s] * nyy[s];

        I1s[s] = (double *) xmalloc(sizes * sizeof(double));
        I2s[s] = (double *) xmalloc(sizes * sizeof(double));
        us[s]  = (double *) xmalloc(sizes * sizeof(double));
        vs[s]  = (double *) xmalloc(sizes * sizeof(double));

        // Compute the zoom from the previous finer scale.
        zoom_out(I1s[s-1], I1s[s], nxx[s-1], nyy[s-1], scaleFactor);
        zoom_out(I2s[s-1], I2s[s], nxx[s-1], nyy[s-1], scaleFactor);
    }

    // Initialize the OF arrays.
    for (i=0; i < nxx[nScales-1] * nyy[nScales-1]; i++) {
        us[nScales-1][i] = 0.0;
        vs[nScales-1][i] = 0.0;
    }

    // Pyramidal approximation to the OF.
    for (s=nScales-1; s>=0; s--) {
        if (verbose)
            printf("Scale: %d %dx%d, nScales=%d\n",
			s, nxx[s], nyy[s], nScales);

        image2warped = (double *) xmalloc(nxx[s] * nyy[s] * sizeof(double));

        // Warp the second image before computing the derivatives.
        bicubic_interpolation_warp(I2s[s], us[s], vs[s], image2warped,
                                   nxx[s], nyy[s], false);

        // Compute the optical flow.
        calcCLG_OF(I1s[s], image2warped, us[s], vs[s], nyy[s], nxx[s],
                   iterations, alpha, rho, wFactor, verbose, coupledMode);

        free(image2warped);

        // If this was the last scale, finish now.
        if (!s) break;
        // Otherwise, upsample the optical flow.

        // Zoom the OF for the next finer scale.
        zoom_in(us[s], us[s-1], nxx[s], nyy[s], nxx[s-1], nyy[s-1]);
        zoom_in(vs[s], vs[s-1], nxx[s], nyy[s], nxx[s-1], nyy[s-1]);

        // Scale the OF with the appropriate zoom factor.
        for (i=0; i < nxx[s-1] * nyy[s-1]; i++) {
            us[s-1][i] *= 1.0 / scaleFactor;
            vs[s-1][i] *= 1.0 / scaleFactor;
        }
    }

    // Free the allocated memory.
    for (i=1; i<nScales; i++) {
        if (verbose)
            printf("Free Scale: %d %dx%d\n", s, nxx[s], nyy[s]);

        free(I1s[i]);
        free(I2s[i]);
        free(us[i]);
        free(vs[i]);
    }

    if (verbose)
        printf("calcMSCLG_OF: done\n");

    return 1;
}

#endif //CLG_OF_C
