// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
// <jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
// All rights reserved.

/**
 * @file clg_of_main.c
 * @brief Entry point of CLG optical flow demo program.
 * @author Jorge Jara <jjara@dcc.uchile.cl>
 */
#include "iio.h"
#include "clg_of.h"


/**
 * read_image
 *
 * Reads an image given its path+name. Returns -1 if the image does not exist.
 *
 * Parameters:
 *
 * fname     name of the input image to read.
 * f         output matrix for the image.
 * w         output width of the image.
 * h         output height of the image.
 * pixeldimg output number of channels of the image.
 *
 */
static int read_image(const char *fname,
                      float **f,
                      int *w,
                      int *h,
                      int *pixeldim) {
    *f = iio_read_image_float_vec(fname, w, h, pixeldim);
    return *f ? 0 : -1;
}


/**
 * main
 *
 * Demo for CLG OF computation, main function. The output is a 2-channel
 * TIFF image with float values, representing the flow field (u,v).
 *
 * The following parameters are expected from the console:
 *
 *  -input_image_1 past image.
 *  -input_image_1 current image.
 *  -alpha         global regularization coefficient
 *  -rho           local smoothing coefficient.
 *  -sigma         standard deviation for gaussian smoothing of the input images
 *                 (optional pre-processing).
 *  -numIt         maximum iteration number.
 *  -w             SOR relaxation factor, between 0 and 2. Default: 1.9.
 *  -nScales       number of scales for the pyramidal approach. Default: 1.
 *  -scaleFactor   scale factor, between 0 and 1. Default: 0.5.
 *  -coupledMode   iteration type: 1->Gauss-Seidel, 0->SOR. Default: 1.
 *  -verbose       shows (1) or hide (0) messages. Default: 1.
 *  -output        output OF field (2-channel TIFF image with float values).
 *
 */
int main(int argc, char *argv[]) {

    if (argc != 8 && argc != 13) {

        fprintf(stderr, "Usage: %s input_image_1 input_image_2 alpha rho "
                "sigma numIt [w nScales zoom_factor coupledMode verbose] "
                "output\n", argv[0]);
    } else {

        int w, h, pixeldim;

        // Input file names.
        const char *image1 = argv[1];
        const char *image2 = argv[2];

        // Parameters for CLG OF calculation.
        double alpha = (double) atof(argv[3]);
        double rho   = (double) atof(argv[4]);
        double sigma = (double) atof(argv[5]);
        int numIt    = atoi(argv[6]);

        // Output file name.
        char *outputFile;

        // Default values for multiscale parameters.
        double wFactor     = 1.9;
        int nScales        = 1;
        double scaleFactor = 0.65;
        int coupledMode    = 1;

        // Verbose mode.
        int verbose = 1;

        if (argc == 8) {

            outputFile = argv[7];

        } else { // 12 input arguments (optical parameters are provided).

            wFactor     = atof(argv[7]);
            nScales     = atoi(argv[8]);
            scaleFactor = (double) atof(argv[9]);
            coupledMode = atoi(argv[10]);
            verbose     = atoi(argv[11]);
            outputFile  = argv[12];
        }

        if (verbose) {
            printf("Parameters: alpha=%f, rho=%f, ", alpha, rho);
            printf("sigma=%f, numIt=%d, w=%f, ", sigma, numIt, wFactor);
            printf("nScales=%d, zoom_factor=%f\n", nScales, scaleFactor);
            printf("Verbose Mode=1, Coupled Mode=%d\n", coupledMode);
            printf("Input: %s, %s, Output: %s\n\n", image1, image2, outputFile);
            printf("allocating memory for data structures\n");
        }

        double* of_u;
        double* of_v;
        double* of_im1;
        double* of_im2;

        float *im1;
        float *im2;

        if (read_image(image1, &im1, &w, &h, &pixeldim) != -1 &&
            read_image(image2, &im2, &w, &h, &pixeldim) != -1) {

           if (verbose)
                printf("Two images loaded:\n\tim1: %dx%d image with %d "
                       "channel(s)\n \tim2, %dx%d image with %d channel(s)\n",
                       w, h, pixeldim, w, h, pixeldim);

            // Memory allocation for variables that require it.
            of_u   = (double *) malloc((size_t) h * w * sizeof(double));
            of_v   = (double *) malloc((size_t) h * w * sizeof(double));
            of_im1 = (double *) malloc((size_t) h * w * sizeof(double));
            of_im2 = (double *) malloc((size_t) h * w * sizeof(double));

            if (verbose)
                printf("copying and casting data for double precision\n");

            int i = 0;
            int j = 0;

            for (i=0; i<h; i++) {
                for (j=0; j<w; j++) {
                    of_im1[w*i+j] = (double) im1[w*i+j];
                    of_im2[w*i+j] = (double) im2[w*i+j];
                    of_u[w*i+j] = 0.0;
                    of_v[w*i+j] = 0.0;
                }
            }

            // CLG OF computation.
            if (verbose)
                printf("call to calc_clg\n");

            // At the highest pyramid level the image must be have
            // enough pixels to apply the smoothing filter.
            double maxSz = DEFAULT_GAUSSIAN_WINDOW_SIZE*rho + 1.0;
            int nScalesMax = (int) floor(1.0 + log(maxSz/fmin(1.0*w, 1.0*h))
                                   / log(scaleFactor));

            if (nScales > nScalesMax) {
                nScales = nScalesMax;
                if (verbose)
                    printf("nScales corrected to the max value allowed: %d\n",
                           nScales);
            }

            int res = calcMSCLG_OF(of_im1, of_im2, of_u, of_v, h, w, numIt,
                                   alpha, rho, sigma, wFactor, nScales,
                                   scaleFactor, coupledMode, verbose);

            // Save results.
            float* x = (float *) malloc((size_t) w * h * 2 * sizeof(float));

            FILE *vxfile;
            FILE *vyfile;

            if (verbose) {
                char filename_vx[255];
                char filename_vy[255];

                if (coupledMode == 0) {
                    sprintf(filename_vx, "vx_sor.txt");
                    sprintf(filename_vy, "vy_sor.txt");
                } else {
                    sprintf(filename_vx, "vx_coupled.txt");
                    sprintf(filename_vy, "vy_coupled.txt");
                }

                printf("Filename vx, %s\n", filename_vx); fflush(stdout);
                printf("Filename vy, %s\n", filename_vy); fflush(stdout);

                vxfile = fopen(filename_vx, "w");
                vyfile = fopen(filename_vy, "w");
            }

            for (i=0; i<h; i++) {
                for (j=0; j<w; j++) {

                    x[      w*i+j] = of_u[w*i+j];
                    x[w*h + w*i+j] = of_v[w*i+j];
                    if (verbose) {
                        fprintf(vxfile, "%f ", of_u[w*i+j]);
                        fprintf(vyfile, "%f ", of_v[w*i+j]);
                    }
                }
                if (verbose) {
                    fprintf(vxfile, "\n");
                    fprintf(vyfile, "\n");
                }
            }

            if (verbose) {
                fclose(vxfile);
                fclose(vyfile);
            }

            iio_save_image_float_split(outputFile, x, w, h, 2);

            // Free memory.
            free(x);
            free(im1);
            free(im2);
            free(of_im1);
            free(of_im2);
            free(of_u);
            free(of_v);

            if (verbose)
                printf("MS CLG optical flow computation done.\n");
        } else {
            fprintf(stderr, "input images cannot be found!.\n");
        }
        return 0;
    }
}
