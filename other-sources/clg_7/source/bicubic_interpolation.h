// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BICUBIC_INTERPOLATION_H
#define BICUBIC_INTERPOLATION_H

#define BOUNDARY_CONDITION 0

#include <stdbool.h>

static int neumann_bc(int x, int nx, bool *out);

static int periodic_bc(int x, int nx, bool *out);

static int symmetric_bc(int x, int nx, bool *out);

static double cubic_interpolation_cell (
    double v[4],  //interpolation points
    double x);    //point to be interpolated

static double bicubic_interpolation_cell (
    double p[4][4], //array containing the interpolation points
    double x,       //x position to be interpolated
    double y);      //y position to be interpolated

double bicubic_interpolation_at(
    const double *input, //image to be interpolated
    const double  uu,    //x component of the vector field
    const double  vv,    //y component of the vector field
    const int    nx,    //image width
    const int    ny,    //image height
    bool         border_out); //if true, return zero outside the region

void bicubic_interpolation_warp(
    const double *input,  //image to be warped
    const double *u,      //x component of the vector field
    const double *v,      //y component of the vector field
    double       *output, //warped output image with bicubic interpolation
    const int    nx,     //image width
    const int    ny,     //image height
    bool         border_out);//if true, put zeros outside the region


#endif//BICUBIC_INTERPOLATION_H
