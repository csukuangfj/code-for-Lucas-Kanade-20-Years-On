// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
//<jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>
// All rights reserved.

#ifndef UTIL_C
#define UTIL_C

/**
 * @file util.c
 * @brief Utility function for the CLG optical flow for 2D grayscale images.
 * @author Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano <jdelpian@uandes.cl>
 * and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
 */
#include "util.h"


/**
 * lin2by2det
 *
 * Determinant for a 2x2 matrix, assuming the form
 * a b
 * c d
 *
 */
double lin2by2det(double a, double b, double c, double d) {
    return a*d - b*c;
}


/**
 * pMatrix
 *
 * Allocates memory for a 2-dimensional array (matrix).
 *
 * Parameters:
 *
 * nRows The number of rows for the requested matrix (first dimension).
 * nCols The number of columns for the requested matrix (second dimension).
 *
 */
double **pMatrix(int nRows, int nCols) {

    double **mat = (double **) malloc((size_t) (nRows * sizeof(double*)));

    if (!mat)
        exit(fprintf(stderr, "out of memory\n"));

    int i = 0;
    for (i=0; i<nRows; i++) {
        mat[i] = (double *) malloc((size_t) (nCols * sizeof(double)));
        if (!mat[i])
            exit(fprintf(stderr, "out of memory\n"));
    }
    return mat;
}


/**
 * freePmatrix
 *
 * Free memory for a given 2-dimensional array (matrix).
 *
 * Parameters:
 *
 * mat   Pointer to the array to be freed.
 * nRows The number of rows for the requested matrix (first dimension).
 *
 */
void freePmatrix(double **mat, int nRows) {

    int i=0, j=0;

    for (i=0; i<nRows; i++)
        free(mat[i]);

    free(mat);
}


/**
 * correctIndex
 *
 * Verifies if a given index is in the range [0,size-1] and corrects it, if
 * necessary, using the closest number from [0,size-1].
 *
 * p    The index to check.
 * size Size of the dimension.
 */
int correctIndex(int p, int size) {
    return (int) fmin(fmax(0.0, (float) p), (float)(size-1.0));
}


/**
 * computeDerivatives
 *
 * Discrete derivative computations for the input time frames (images). Spatial
 * derivatives are computed using the second image as the "current" time frame.
 *
 * Parameters:
 *
 * image1  Pointer to the first time frame (image).
 * image2  Pointer to the second time frame (image).
 * dfdx    Pointer to the x-coordinate spatial derivatives matrix.
 * dfdy    Pointer to the y-coordinate spatial derivatives matrix.
 * nRows   Number of rows of the images.
 * nCols   Number of columns of the images.
 * verbose Sets the verbose mode on or off.
 *
 */
void computeDerivatives(double **image1,
                        double **image2,
                        double **dfdx,
                        double **dfdy,
                        int nRows,
                        int nCols,
                        int verbose) {
    if (verbose)
        printf("  computeDerivatives\n");

    int i=1, j=1;

    // Using a 3x1 discrete derivative stencil according to Bruhn et al.
    // The operator derives from a 2nd order Taylor expansion.
    double imDerStencil[] = {-0.5, 0.0, 0.5};

    int imDerStencilSize = 3;
    int imDerStencilSide = imDerStencilSize / 2;

    // The spatial gradient can be computed in image1 or image2.
    double **gradientIm = image1;

    if (verbose)
        printf("    loop for derivatives computation: kernel size %d"
               " (side %d)\n", imDerStencilSize, imDerStencilSide);

    for (i=0; i<nRows; i++) {

        for (j=0; j<nCols; j++) {

            dfdx[i][j] = 0.0;
            dfdy[i][j] = 0.0;

            int count = 0;
            for (count=0; count<imDerStencilSize; count++) {
                 dfdy[i][j] += imDerStencil[count]*
                   gradientIm[correctIndex(i+count-imDerStencilSide, nRows)][j];
            }
            count = 0;
            for (count=0; count<imDerStencilSize; count++) {
                 dfdx[i][j] += imDerStencil[count]*
                   gradientIm[i][correctIndex(j+count-imDerStencilSide, nCols)];
            }

        }
    }

    if (verbose)
        printf("    done.\n");
}


/**
 * computeJTensor
 *
 * Computes and sets the J tensor values.
 *
 * Parameters:
 *
 * dfdx            Pointer to the x-coordinate spatial derivatives input matrix.
 * dfdy            Pointer to the y-coordinate spatial derivatives input matrix.
 * dfdt            Pointer to the time derivatives input matrix.
 * J[JROWS][JCOLS] Pointer to the motion tensor matrix J, computed at each image
 *                 pixel, containing its corresponding discrete derivatives.
 *                 For a given pixel [i,j], the tensor has the form
 *
 *                 J[0][0][i][j] = dfdx[i][j] * dfdx[i][j];
 *                 J[0][1][i][j] = dfdx[i][j] * dfdy[i][j];
 *                 J[0][2][i][j] = dfdx[i][j] * dfdt[i][j];
 *                 J[1][1][i][j] = dfdy[i][j] * dfdy[i][j];
 *                 J[1][2][i][j] = dfdy[i][j] * dfdt[i][j];
 *                 J[2][2][i][j] = dfdt[i][j] * dfdt[i][j];
 *
 *                 with f being the second (current) frame.
 *                 Since the tensor is symmetric, i.e.
 *
 *                 J[1][0][i][j] = J[0][1][i][j];
 *                 J[2][0][i][j] = J[0][2][i][j];
 *                 J[2][1][i][j] = J[1][2][i][j];
 *
 *                 only the upper half and the diagonal are stored.
 * nRows           Number of rows of the input derivative matrices.
 * nCols           Number of columns of the input derivative matrices.
 *
 */
void computeJTensor(double **dfdx,
                    double **dfdy,
                    double **dfdt,
                    double **J[JROWS][JCOLS],
                    int nRows,
                    int nCols) {
    int i, j;

    for (i=0; i<nRows; i++)
        for (j=0; j<nCols; j++) {
            J[0][0][i][j] = dfdx[i][j] * dfdx[i][j];
            J[0][1][i][j] = dfdx[i][j] * dfdy[i][j];
            J[0][2][i][j] = dfdx[i][j] * dfdt[i][j];
            J[1][1][i][j] = dfdy[i][j] * dfdy[i][j];
            J[1][2][i][j] = dfdy[i][j] * dfdt[i][j];
            J[2][2][i][j] = dfdt[i][j] * dfdt[i][j];
        }
}


/**
 * matrixSmooth
 *
 * Performs a gaussian smoothing on a given matrix, with a specified size for
 * the kernel (filter) and standard deviation value.
 *
 * Parameters:
 * matrix      Pointer to the input matrix.
 * nRows       Number of rows of the input matrix.
 * nCols       Number of columns of the input matrix.
 * kernelSigma Standard deviation for the gaussian kernel.
 *
 */
void matrixSmooth(double **matrix,
                  int nRows,
                  int nCols,
                  double kernelSigma) {
    int i, j;

    // Convert from matrix to array form.
    double *matrixArray = xmalloc(nRows * nCols * sizeof(double));

    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            matrixArray[j + i*nCols] = matrix[i][j];
        }
    }

    // Convolve.
    gaussian(matrixArray, nCols, nRows, kernelSigma);

    // Reassign result values to the matrix.
    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            matrix[i][j] = matrixArray[j + i*nCols];
        }
    }

    // free memory
    free(matrixArray);
}


/**
 * xmalloc
 *
 * This function is like "malloc", but always returns a valid pointer.
 *
 * Parameters:
 * size The amount of memory to allocate, passed to the "malloc" function.
 *
 */
void *xmalloc(size_t size) {

    void *p = malloc(size);

    if (!p)
        exit(fprintf(stderr, "out of memory\n"));

    return p;
}

#endif
