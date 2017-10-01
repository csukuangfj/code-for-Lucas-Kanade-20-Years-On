// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_OPTIC_FLOW_H
#define BROX_OPTIC_FLOW_H

#include <vector>
#include <iostream>

#include "mask.h"
#include "zoom.h"
#include "bicubic_interpolation.h"

#define EPSILON 0.001
#define MAXITER 300
#define SOR_PARAMETER 1.9
#define GAUSSIAN_SIGMA 0.8

/**
  *
  * Compute the coefficients of the robust functional (data term)
  *
**/
void psi_data(
    const float *I1,  //first image
    const float *I2,  //second image
    const float *I2x, //gradient of the second image
    const float *I2y, //gradient of the second image
    const float *du,  //motion increment
    const float *dv,  //motion increment
    float *psip,      //output coefficients
    const int nx,     //image width
    const int ny      //image height
)
{
    const int size = nx * ny;

    //compute 1/(sqrt((I2-I1+I2x*du+I2y*dv)²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float dI  = I2[i] - I1[i] + I2x[i] * du[i] + I2y[i] * dv[i];
	const float dI2 = dI * dI;

	psip[i] = 1. / sqrt(dI2 + EPSILON * EPSILON);
    }
}


/**
  *
  * Compute the coefficients of the robust functional (gradient term)
  *
**/
void psi_gradient(
    const float *I1x,  //gradient of the first image
    const float *I1y,  //gradient of the first image
    const float *I2x,  //gradient of the second image
    const float *I2y,  //gradient of the second image
    const float *I2xx, //second derivatives of the second image
    const float *I2xy, //second derivatives of the second image
    const float *I2yy, //second derivatives of the second image
    const float *du,   //motion increment
    const float *dv,   //motion increment
    float *psip,       //output coefficients
    const int nx,      //image width
    const int ny       //image height
)
{
    const int size = nx * ny;

    //compute 1/(sqrt(|DI2-DI1+HI2*(du,dv)|²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float dIx = I2x[i] - I1x[i] + I2xx[i] * du[i] + I2xy[i] * dv[i];
	const float dIy = I2y[i] - I1y[i] + I2xy[i] * du[i] + I2yy[i] * dv[i];
	const float dI2 = dIx * dIx + dIy * dIy;

	psip[i] = 1. / sqrt(dI2 + EPSILON * EPSILON);
    }
}


/**
  *
  * Compute the coefficients of the robust functional (smoothness term)
  *
**/
void psi_smooth(
    const float *ux, //gradient of x component of the optical flow
    const float *uy, //gradient of x component of the optical flow
    const float *vx, //gradient of y component of the optical flow
    const float *vy, //gradient of y component of the optical flow
    float *psi,      //output coefficients
    const int nx,    //image width
    const int ny     //image height
)
{
    const int size = nx * ny;

    //compute 1/(sqrt(ux²+uy²+vx²+vy²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float du  = ux[i] * ux[i] + uy[i] * uy[i];
	const float dv  = vx[i] * vx[i] + vy[i] * vy[i];
	const float d2  = du + dv;

	psi[i] = 1. / sqrt(d2 + EPSILON * EPSILON);
    }
}


/**
 * 
 *  SOR iteration in one position
 * 
 */
inline float sor_iteration(
    const float *Au,   //constant part of the numerator of u
    const float *Av,   //constant part of the numerator of v
    const float *Du,   //denominator of u
    const float *Dv,   //denominator of v
    const float *D,    //constant part of the numerator
    float       *du,   //x component of the motion increment
    float       *dv,   //y component of the motion increment
    const float alpha, //alpha smoothness parameter
    const float *psi1, //coefficients of the divergence
    const float *psi2, 
    const float *psi3,
    const float *psi4,
    const int   i,     //current row
    const int   i0,    //previous row
    const int   i1,    //following row
    const int   j,     //current column
    const int   nx,    //number of columns
    const int   j0,    //previous column
    const int   j1     //following column
)
{
    //set the SOR extrapolation parameter
    const float w = SOR_PARAMETER;

    //calculate the position in the array
    const int k = i * nx + j;

    //compute the divergence part of the numerator
    const float div_du = psi1[k] * du[k+i1] + psi2[k] * du[k-i0] + 
			  psi3[k] * du[k+j1] + psi4[k] * du[k-j0] ;
    const float div_dv = psi1[k] * dv[k+i1] + psi2[k] * dv[k-i0] + 
			  psi3[k] * dv[k+j1] + psi4[k] * dv[k-j0] ;

    const float duk = du[k];
    const float dvk = dv[k];

    //update the motion increment
    du[k] = (1.-w) * du[k] + w * (Au[k] - D[k] * dv[k] + alpha * div_du) / Du[k];
    dv[k] = (1.-w) * dv[k] + w * (Av[k] - D[k] * du[k] + alpha * div_dv) / Dv[k];

    //return the covergence error in this position
    return (du[k] - duk) * (du[k] - duk) + (dv[k] - dvk) * (dv[k] - dvk); 
}


/**
  *
  * Compute the optic flow with the Brox spatial method
  *
**/
void brox_optic_flow
(
    const float *I1,         //first image
    const float *I2,         //second image
    float *u, 		      //x component of the optical flow
    float *v, 		      //y component of the optical flow
    const int    nx,         //image width
    const int    ny,         //image height
    const float  alpha,      //smoothness parameter
    const float  gamma,      //gradient term parameter
    const float  TOL,        //stopping criterion threshold
    const int    inner_iter, //number of inner iterations
    const int    outer_iter, //number of outer iterations
    const bool   verbose     //switch on messages
)
{
    const int size = nx * ny;

    //allocate memory
    float *du    = new float[size];
    float *dv    = new float[size];

    float *ux    = new float[size];
    float *uy    = new float[size];
    float *vx    = new float[size];
    float *vy    = new float[size];

    float *I1x   = new float[size];
    float *I1y   = new float[size];
    float *I2x   = new float[size];
    float *I2y   = new float[size];
    float *I2w   = new float[size];
    float *I2wx  = new float[size];
    float *I2wy  = new float[size];
    float *I2xx  = new float[size];
    float *I2yy  = new float[size];
    float *I2xy  = new float[size];
    float *I2wxx = new float[size];
    float *I2wyy = new float[size];
    float *I2wxy = new float[size];

    float *div_u = new float[size];
    float *div_v = new float[size];
    float *div_d = new float[size];

    float *Au    = new float[size];
    float *Av    = new float[size];
    float *Du    = new float[size];
    float *Dv    = new float[size];
    float *D     = new float[size];

    float *psid  = new float[size];
    float *psig  = new float[size];
    float *psis  = new float[size];
    float *psi1  = new float[size];
    float *psi2  = new float[size];
    float *psi3  = new float[size];
    float *psi4  = new float[size];

    //compute the gradient of the images
    gradient(I1, I1x, I1y, nx, ny);
    gradient(I2, I2x, I2y, nx, ny);

    //compute second order derivatives
    Dxx(I2, I2xx, nx, ny);
    Dyy(I2, I2yy, nx, ny);
    Dxy(I2, I2xy, nx, ny);

    //outer iterations loop
    for(int no = 0; no < outer_iter; no++)
    {
	//warp the second image and its derivatives
	bicubic_interpolation(I2,   u, v, I2w,   nx, ny, true);
	bicubic_interpolation(I2x,  u, v, I2wx,  nx, ny, true);
	bicubic_interpolation(I2y,  u, v, I2wy,  nx, ny, true);
	bicubic_interpolation(I2xx, u, v, I2wxx, nx, ny, true);
	bicubic_interpolation(I2xy, u, v, I2wxy, nx, ny, true);
	bicubic_interpolation(I2yy, u, v, I2wyy, nx, ny, true);

	//compute the flow gradient
	gradient(u, ux, uy, nx, ny);
	gradient(v, vx, vy, nx, ny);

	//compute robust function Phi for the smoothness term
	psi_smooth(ux, uy, vx, vy, psis, nx, ny);

	//compute coefficients of Phi functions in divergence
	psi_divergence(psis, psi1, psi2, psi3, psi4, nx, ny);

	//compute the divergence for the gradient of w
	divergence_u(u, v, psi1, psi2, psi3, psi4, div_u, div_v, nx, ny);

	for(int i = 0; i < size; i++)
	{
	    //compute the coefficents of dw[i] in the smoothness term
	    div_d[i] = alpha * (psi1[i] + psi2[i] + psi3[i] + psi4[i]);

	    //initialize the motion increment
	    du[i] = dv[i] = 0;
	}

	//inner iterations loop
	for(int ni = 0; ni < inner_iter; ni++)
	{
	    //compute robust function Phi for the data and gradient terms
	    psi_data(I1, I2w, I2wx, I2wy, du, dv,  psid, nx, ny);
	    psi_gradient(I1x, I1y, I2wx, I2wy, I2wxx, I2wxy, I2wyy, du, dv, psig, nx, ny);

	    //store constant parts of the numerical scheme
	    for(int i = 0; i < size; i++)
	    {
		const float p = psid[i];
		const float g = gamma * psig[i];

		//brightness constancy term
		const float dif = I2w[i] - I1[i];
		const float BNu = -p * dif * I2wx[i];
		const float BNv = -p * dif * I2wy[i];
		const float BDu = p * I2wx[i] * I2wx[i];
		const float BDv = p * I2wy[i] * I2wy[i];

		//gradient constancy term
		const float dx  = (I2wx[i] - I1x[i]);
		const float dy  = (I2wy[i] - I1y[i]);
		const float GNu = -g * (dx * I2wxx[i] + dy * I2wxy[i]);
		const float GNv = -g * (dx * I2wxy[i] + dy * I2wyy[i]);
		const float GDu =  g * (I2wxx[i] * I2wxx[i] + I2wxy[i] * I2wxy[i]);
		const float GDv =  g * (I2wyy[i] * I2wyy[i] + I2wxy[i] * I2wxy[i]);
		const float DI  = (I2wxx[i] + I2wyy[i]) * I2wxy[i];
		const float Duv =  p * I2wy[i] * I2wx[i] + g * DI;	    

		Au[i] = BNu + GNu + alpha * div_u[i];
		Av[i] = BNv + GNv + alpha * div_v[i];
		Du[i] = BDu + GDu + div_d[i];
		Dv[i] = BDv + GDv + div_d[i];
		D [i] = Duv;
	    }

	    //sor iterations loop
	    float error = 1000;
	    int nsor = 0;
	    
	    while( error > TOL && nsor < MAXITER)
	    {
		error = 0;
		nsor++;
		
		//update the motion increment in the center of the images
		for(int i = 1; i < ny-1; i++)
		    for(int j = 1; j < nx-1; j++) {

			error += sor_iteration(
			      Au, Av, Du, Dv, D, du, dv, alpha,  
			      psi1, psi2, psi3, psi4,
			      i, nx, nx, j, nx, 1, 1
			);
		    }

		//update the motion increment in the first and last rows
		for(int j = 1; j < nx-1; j++)
		{
		    error += sor_iteration(
			Au, Av, Du, Dv, D, du, dv, alpha,  
			psi1, psi2, psi3, psi4,
			0, 0, nx, j, nx, 1, 1
		    );    

		    error += sor_iteration(
			Au, Av, Du, Dv, D, du, dv, alpha,  
			psi1, psi2, psi3, psi4,
			ny-1, nx, 0, j, nx, 1, 1
		    );    
		}

		//update the motion increment in the first and last columns
		for(int i = 1; i < ny-1; i++)
		{
		    error += sor_iteration(
			Au, Av, Du, Dv, D, du, dv, alpha,  
			psi1, psi2, psi3, psi4,
			i, nx, nx, 0, nx, 0, 1
		    );

		    error += sor_iteration(
			Au, Av, Du, Dv, D, du, dv, alpha,  
			psi1, psi2, psi3, psi4,
			i, nx, nx, nx-1, nx, 1, 0
		    );
		}
		
		//process the top-left corner (0,0)
		error += sor_iteration(
		    Au, Av, Du, Dv, D, du, dv, alpha,  
		    psi1, psi2, psi3, psi4,
		    0, 0, nx, 0, nx, 0, 1
		);

		//process the top-right corner (0,nx-1)
		error += sor_iteration(
		    Au, Av, Du, Dv, D, du, dv, alpha,  
		    psi1, psi2, psi3, psi4,
		    0, 0, nx, nx-1, nx, 1, 0
		);

		//process the bottom-left corner (ny-1,0)
		error += sor_iteration(
		    Au, Av, Du, Dv, D, du, dv, alpha,  
		    psi1, psi2, psi3, psi4,
		    ny-1, nx, 0, 0, nx, 0, 1
		);

		//process the bottom-right corner (ny-1,nx-1)
		error += sor_iteration(
		    Au, Av, Du, Dv, D, du, dv, alpha,  
		    psi1, psi2, psi3, psi4,
		    ny-1, nx, 0, nx-1, nx, 1, 0
		);
		
		error = sqrt(error / size);
	    }
	    
	    if(verbose) std::cout << "Iterations: " << nsor << std::endl; 
	}

	//update the flow with the estimated motion increment
	for(int i = 0; i < size; i++)
	{
	    u[i] += du[i];
	    v[i] += dv[i];
	}
    }

    //delete allocated memory
    delete []du;
    delete []dv;

    delete []ux;
    delete []uy;
    delete []vx;
    delete []vy;

    delete []I1x;
    delete []I1y;
    delete []I2x;
    delete []I2y;
    delete []I2w;
    delete []I2wx;
    delete []I2wy;
    delete []I2xx;
    delete []I2yy;
    delete []I2xy;
    delete []I2wxx;
    delete []I2wyy;
    delete []I2wxy;

    delete []div_u;
    delete []div_v;    
    delete []div_d;

    delete []Au;
    delete []Av;
    delete []Du;
    delete []Dv;
    delete []D;

    delete []psid;
    delete []psig;
    delete []psis;
    delete []psi1;
    delete []psi2;
    delete []psi3;
    delete []psi4;
}


/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
void image_normalization(
    const float *I1,   //input image 1
    const float *I2,   //input image 2
    float       *I1n,  //normalized output image 1
    float       *I2n,  //normalized output image 2 
    int          size  //size of the image
)
{
    //compute the max and min values of the images
    const float max0 = *std::max_element(I1, &I1[size]);
    const float max1 = *std::max_element(I2, &I2[size]);
    const float min0 = *std::min_element(I1, &I1[size]);
    const float min1 = *std::min_element(I2, &I2[size]);

    //compute the global max and min
    const float max = std::max(max0, max1);
    const float min = std::min(min0, min1);
    const float den = max - min;

    if(den > 0)
	//normalize the images between 0 and 255
	for(int i = 0; i < size; i++)
	{
	    I1n[i] = 255.0 * (I1[i] - min) / den;
	    I2n[i] = 255.0 * (I2[i] - min) / den;
	}

    else
	//copy the original data
	for(int i = 0; i < size; i++)
	{
	    I1n[i] = I1[i];
	    I2n[i] = I2[i];
	}
}


/**
  *
  *  Multiscale approach for computing the optical flow
  *
**/
void brox_optic_flow(
    const float *I1,         //first image
    const float *I2,         //second image
    float *u, 		      //x component of the optical flow
    float *v, 		      //y component of the optical flow
    const int    nxx,        //image width
    const int    nyy,        //image height
    const float  alpha,      //smoothness parameter
    const float  gamma,      //gradient term parameter
    const int    nscales,    //number of scales
    const float  nu,         //downsampling factor
    const float  TOL,        //stopping criterion threshold
    const int    inner_iter, //number of inner iterations
    const int    outer_iter, //number of outer iterations
    const bool   verbose     //switch on messages
)
{
    int size = nxx * nyy;

    std::vector<float *> I1s(nscales);
    std::vector<float *> I2s(nscales);
    std::vector<float *> us (nscales);
    std::vector<float *> vs (nscales);

    std::vector<int> nx(nscales);
    std::vector<int> ny(nscales);

    I1s[0] = new float[size];
    I2s[0] = new float[size];

    //normalize the input images between 0 and 255
    image_normalization(I1, I2, I1s[0], I2s[0], size);

    //presmoothing the finest scale images
    gaussian(I1s[0], nxx, nyy, GAUSSIAN_SIGMA);
    gaussian(I2s[0], nxx, nyy, GAUSSIAN_SIGMA);

    us [0] = u;
    vs [0] = v;
    nx [0] = nxx;
    ny [0] = nyy;

    //create the scales
    for(int s = 1; s < nscales; s++)
    {
	zoom_size(nx[s-1], ny[s-1], nx[s], ny[s], nu);
	const int sizes = nx[s] * ny[s];

	I1s[s] = new float[sizes];
	I2s[s] = new float[sizes];
	us[s]  = new float[sizes];
	vs[s]  = new float[sizes];

	//compute the zoom from the previous scale
	zoom_out(I1s[s-1], I1s[s], nx[s-1], ny[s-1], nu);
	zoom_out(I2s[s-1], I2s[s], nx[s-1], ny[s-1], nu);
    }

    //initialization of the optical flow at the coarsest scale
    for(int i = 0; i < nx[nscales-1] * ny[nscales-1]; i++)
	us[nscales-1][i] = vs[nscales-1][i] = 0.0;


    //pyramidal approach for computing the optical flow
    for(int s = nscales-1; s >= 0; s--)
    {
	if(verbose) std::cout << "Scale: " << s << std::endl;

	//compute the optical flow for the current scale
	brox_optic_flow(
	    I1s[s], I2s[s], us[s], vs[s], nx[s], ny[s], 
	    alpha, gamma, TOL, inner_iter, outer_iter, verbose
	);

	//if it is not the finer scale, then upsample the optical flow and adapt it conveniently
	if(s)
	{
	    zoom_in(us[s], us[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
	    zoom_in(vs[s], vs[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);

	    for(int i = 0; i < nx[s-1] * ny[s-1]; i++)
	    {
		us[s-1][i] *= 1.0 / nu;
		vs[s-1][i] *= 1.0 / nu;
	    }
	}
    }

    //delete allocated memory
    delete []I1s[0];
    delete []I2s[0];

    for(int i = 1; i < nscales; i++)
    {
	delete []I1s[i];
	delete []I2s[i];
	delete []us [i];
	delete []vs [i];
    }
}

#endif
