// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Jorge Jara <jjara@dcc.uchile.cl>, Jose Delpiano
// <jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>.
// All rights reserved.

/**
 * @file test_rgb2gray.c
 * @brief Utility function for RGB-grayscale image conversion.
 * @author Jorge Jara <jjara@dcc.uchile.cl>
 */
#include "iio.h"
#include <stdio.h>  // Used for "printf".
#include <stdlib.h> // Used for "free".

int main(int argc, char *argv[]) {

    if (argc != 3) {

        printf("Usage: %s input_image output_image\n", argv[0]);

    } else {

        int w, h, pixeldim;
        float *x = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
        printf("Got a %dx%d image with %d channels\n", w, h, pixeldim);

        // conversion to gray scale
        float *xgray = malloc(w * h * sizeof(float));
        if (pixeldim == 3) {
            int i;
            for (i=0; i<w*h; i++) {
                xgray[i] = (6968 * x[pixeldim*i]
                            + 23434 * x[pixeldim*i + 1]
                            +  2366 * x[pixeldim*i + 2]) / 32768;
            }
            // Y = (6968 R + 23434 G + 2366 B) / 32768
        } else {
            xgray = x;
        }

        printf("w = %d, h = %d, pixeldim = %d.\n", w, h, pixeldim);
        iio_save_image_float_vec(argv[2], xgray, w, h, 1);
        free(x);
        return 0;
    }
}
