// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef BROX_OPTIC_FLOW_H
#define BROX_OPTIC_FLOW_H

#include "mask.h"
#include "zoom.h"
#include "bicubic_interpolation.h"

#include <vector>
#include <iostream>
#include <algorithm>

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
    const int size   //images size
)
{
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
    const int size     //images size
)
{
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
    const float *ux, //x derivative of the u component
    const float *uy, //y derivative of the u component
    const float *ut, //t derivative of the u component
    const float *vx, //x derivative of the v component
    const float *vy, //y derivative of the v component
    const float *vt, //t derivative of the v component
    float *psi,      //output coefficients
    const int size   //image size
)
{
    //compute 1/(sqrt(ux²+uy²+vx²+vy²+e²) in each pixel
    for(int i = 0; i < size; i++)
    {
	const float du  = ux[i] * ux[i] + uy[i] * uy[i] + ut[i] * ut[i];
	const float dv  = vx[i] * vx[i] + vy[i] * vy[i] + vt[i] * vt[i] ;
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
    const float *psi5,
    const float *psi6,
    const int   f,     //current frame
    const int   df0,   //previous frame
    const int   df1,   //following frame
    const int   i,     //current row
    const int   ny,    //number of rows
    const int   dy0,   //previous row
    const int   dy1,   //following row
    const int   j,     //current column
    const int   nx,    //number of columns
    const int   dx0,   //previous column
    const int   dx1    //following column
)
{
    //set the SOR extrapolation parameter
    const float w = SOR_PARAMETER;

    //calculate the position in the array
    const int k = f * ny * nx + i * nx + j;

    //compute the divergence part of the numerator
    const float div_du = psi1[k] * du[k+dy1] + psi2[k] * du[k-dy0] + 
			  psi3[k] * du[k+dx1] + psi4[k] * du[k-dx0] +
			  psi5[k] * du[k-df0] + psi6[k] * du[k+df1] ;
    const float div_dv = psi1[k] * dv[k+dy1] + psi2[k] * dv[k-dy0] + 
			  psi3[k] * dv[k+dx1] + psi4[k] * dv[k-dx0] +
			  psi5[k] * dv[k-df0] + psi6[k] * dv[k+df1] ;
      
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
 *  Procedure to compute the motion increment in one frame
 * 
 */
float process_frame(
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
    const float *psi5,
    const float *psi6,
    const int   f,     //current frame
    const int   nx,    //number of columns
    const int   ny,    //number of rows
    const int   df0,   //previous frame
    const int   df1    //following frame   
)
{
    float error = 0;

    //update the motion increment in the center of the images
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	
	    error += sor_iteration(
		Au, Av, Du, Dv, D, du, dv, alpha,  
		psi1, psi2, psi3, psi4, psi5, psi6,
		f, df0, df1, i, ny, nx, nx, j, nx, 1, 1
	    );    
    }

    //update the motion increment in the first and last rows
    for(int j = 1; j < nx-1; j++)
    {
	error += sor_iteration(
	    Au, Av, Du, Dv, D, du, dv, alpha,  
	    psi1, psi2, psi3, psi4, psi5, psi6,
	    f, df0, df1, 0, ny, 0, nx, j, nx, 1, 1
	);    

	error += sor_iteration(
	    Au, Av, Du, Dv, D, du, dv, alpha,  
	    psi1, psi2, psi3, psi4, psi5, psi6,
	    f, df0, df1, ny-1, ny, nx, 0, j, nx, 1, 1
	);    
    }

    //update the motion increment in the first and last columns
    for(int i = 1; i < ny-1; i++)
    {
	error += sor_iteration(
	    Au, Av, Du, Dv, D, du, dv, alpha,  
	    psi1, psi2, psi3, psi4, psi5, psi6,
	    f, df0, df1, i, ny, nx, nx, 0, nx, 0, 1
	);

	error += sor_iteration(
	    Au, Av, Du, Dv, D, du, dv, alpha,  
	    psi1, psi2, psi3, psi4, psi5, psi6,
	    f, df0, df1, i, ny, nx, nx, nx-1, nx, 1, 0
	);
    }


    //process the top-left corner (0,0)
    error += sor_iteration(
	Au, Av, Du, Dv, D, du, dv, alpha,  
	psi1, psi2, psi3, psi4, psi5, psi6,
	f, df0, df1, 0, ny, 0, nx, 0, nx, 0, 1
    );

    //process the top-right corner (0,nx-1)
    error += sor_iteration(
	Au, Av, Du, Dv, D, du, dv, alpha,  
	psi1, psi2, psi3, psi4, psi5, psi6,
	f, df0, df1, 0, ny, 0, nx, nx-1, nx, 1, 0
    );

    //process the bottom-left corner (ny-1,0)
    error += sor_iteration(
	Au, Av, Du, Dv, D, du, dv, alpha,  
	psi1, psi2, psi3, psi4, psi5, psi6,
	f, df0, df1, ny-1, ny, nx, 0, 0, nx, 0, 1
    );

    //process the bottom-right corner (ny-1,nx-1)
    error += sor_iteration(
	Au, Av, Du, Dv, D, du, dv, alpha,  
	psi1, psi2, psi3, psi4, psi5, psi6,
	f, df0, df1, ny-1, ny, nx, 0, nx-1, nx, 1, 0
    );

    return error;
}


/**
  *
  * Compute the optic flow with the Brox temporal method 
  *
**/
void brox_optic_flow
(
    const float *I,          //sequence of images
    float *u, 		      //x component of the optical flow
    float *v, 		      //y component of the optical flow
    const int    nx,         //image width
    const int    ny,         //image height
    const int    frames,     //number of frames
    const float  alpha,      //smoothness parameter
    const float  gamma,      //gradient term parameter
    const float  TOL,        //stopping criterion threshold
    const int    inner_iter, //number of inner iterations
    const int    outer_iter, //number of outer iterations
    const bool   verbose     //switch on messages
)
{
    const int nz    = frames-1;
    const int df    = nx * ny;
    const int size  = df * frames;
    const int size1 = df * nz;

    //allocate memory
    float *du    = new float[size1];
    float *dv    = new float[size1];

    float *ux    = new float[size1];
    float *uy    = new float[size1];
    float *ut    = new float[size1];
    float *vx    = new float[size1];
    float *vy    = new float[size1];
    float *vt    = new float[size1];

    float *Ix   = new float[size];
    float *Iy   = new float[size];
    float *Iw   = new float[size1];
    float *Iwx  = new float[size1];
    float *Iwy  = new float[size1];
    float *Ixx  = new float[size1];
    float *Iyy  = new float[size1];
    float *Ixy  = new float[size1];
    float *Iwxx = new float[size1];
    float *Iwyy = new float[size1];
    float *Iwxy = new float[size1];

    float *div_u = new float[size1];
    float *div_v = new float[size1];
    float *div_d = new float[size1];

    float *Au    = new float[size1];
    float *Av    = new float[size1];
    float *Du    = new float[size1];
    float *Dv    = new float[size1];
    float *D     = new float[size1];

    float *psid = new float[size1];
    float *psig = new float[size1];
    float *psis = new float[size1];
    float *psi1 = new float[size1];
    float *psi2 = new float[size1];
    float *psi3 = new float[size1];
    float *psi4 = new float[size1];
    float *psi5 = new float[size1];
    float *psi6 = new float[size1];

    //compute the gradient of the images
    for(int f = 0; f < frames; f++)
	gradient(&I[f*df], &Ix[f*df], &Iy[f*df], nx, ny);

    //compute second order derivatives
    for(int f = 0; f < nz; f++)
    {
	Dxx(&I[df*(f+1)], &Ixx[f*df], nx, ny);
	Dyy(&I[df*(f+1)], &Iyy[f*df], nx, ny);
	Dxy(&I[df*(f+1)], &Ixy[f*df], nx, ny);
    }

    //outer iterations loop
    for(int no = 0; no < outer_iter; no++)
    {
	//warp the images and their derivatives
	for(int f = 0; f < nz; f++)
	{
	    bicubic_interpolation(&I[df*(f+1)],  &u[f*df], &v[f*df], &Iw[f*df],   nx, ny, true);
	    bicubic_interpolation(&Ix[df*(f+1)], &u[f*df], &v[f*df], &Iwx[f*df],  nx, ny, true);
	    bicubic_interpolation(&Iy[df*(f+1)], &u[f*df], &v[f*df], &Iwy[f*df],  nx, ny, true);
	    bicubic_interpolation(&Ixx[f*df],    &u[f*df], &v[f*df], &Iwxx[f*df], nx, ny, true);
	    bicubic_interpolation(&Ixy[f*df],    &u[f*df], &v[f*df], &Iwxy[f*df], nx, ny, true);
	    bicubic_interpolation(&Iyy[f*df],    &u[f*df], &v[f*df], &Iwyy[f*df], nx, ny, true);
	}
	
	//compute the flow gradient
	gradient3(u, ux, uy, ut, nx, ny, nz);
	gradient3(v, vx, vy, vt, nx, ny, nz);

	//compute robust function Phi for the brightness and gradient terms
	psi_smooth(ux, uy, ut, vx, vy, vt, psis, size1);

	//copmute coefficients of Phi functions in divergence
	psi_divergence(psis, psi1, psi2, psi3, psi4, psi5, psi6, nx, ny, nz);

	//compute the divergence for the gradient of w
	divergence_u(u, v, psi1, psi2, psi3, psi4, psi5, psi6, div_u, div_v, nx, ny, nz);


	for(int i = 0; i < size1; i++)
	{
	    //compute the coefficents of dw[i] in the smoothness term
	    div_d[i] = alpha * (psi1[i] + psi2[i] + psi3[i] + psi4[i] + psi5[i] + psi6[i]);

	    //initialize the motion increment
	    du[i] = dv[i] = 0;
	}


	//inner iterations loop
	for(int ni = 0; ni < inner_iter; ni++)
	{
	    //compute robust function Phi for the brightness and gradient terms
	    psi_data(I, Iw, Iwx, Iwy, du, dv, psid, size1);
	    psi_gradient(Ix, Iy, Iwx, Iwy, Iwxx, Iwxy, Iwyy, du, dv, psig, size1);

	    //store constant parts of the numerical scheme
	    for(int i = 0; i < size1; i++)
	    {
		const float p = psid[i];
		const float g = gamma * psig[i];

		//brightness constancy term
		const float dif = Iw[i] - I[i];
		const float BNu = -p * dif * Iwx[i];
		const float BNv = -p * dif * Iwy[i];
		const float BDu =  p * Iwx[i] * Iwx[i];
		const float BDv =  p * Iwy[i] * Iwy[i];

		//gradient constancy term
		const float dx  = (Iwx[i] - Ix[i]);
		const float dy  = (Iwy[i] - Iy[i]);
		const float GNu = -g * (dx * Iwxx[i] + dy * Iwxy[i]);
		const float GNv = -g * (dx * Iwxy[i] + dy * Iwyy[i]);
		const float GDu =  g * (Iwxx[i] * Iwxx[i] + Iwxy[i] * Iwxy[i]);
		const float GDv =  g * (Iwyy[i] * Iwyy[i] + Iwxy[i] * Iwxy[i]);
		const float DI  = (Iwxx[i] + Iwyy[i]) * Iwxy[i];
		const float Duv =  p * Iwy[i] * Iwx[i] + g * DI;

		Au[i] = BNu + GNu + alpha * div_u[i];
		Av[i] = BNv + GNv + alpha * div_v[i];
		Du[i] = BDu + GDu + div_d[i];
		Dv[i] = BDv + GDv + div_d[i];
		D [i] = Duv;
	    }


	    //sor iterations loop
	    float error = 1000;
	    int   nsor  = 0;	    
	    while(error > TOL && nsor < MAXITER)
	    {
		error = 0;
		nsor++;

		//compute the interior optical flows
		for(int f = 1; f < nz-1; f++)
		{
		    error += process_frame(
			Au, Av, Du, Dv, D, du, dv, alpha, 
			psi1, psi2, psi3, psi4, psi5, psi6,
			f, nx, ny, df, df
		    );    
		} 

		//compute the first optical flow
		error += process_frame(
		    Au, Av, Du, Dv, D, du, dv, alpha, 
		    psi1, psi2, psi3, psi4, psi5, psi6,
		    0, nx, ny, 0, df
		);    

		//compute the last optical flow
		error += process_frame(
		    Au, Av, Du, Dv, D, du, dv, alpha,  
		    psi1, psi2, psi3, psi4, psi5, psi6,
		    nz-1, nx, ny, df, 0
		);
		
		error = sqrt(error / size1);
	    }
	
	    if(verbose) std::cout << "Iterations: " << nsor << std::endl; 
	}
	
	//update the flow with the estimated motion increment
	for(int i = 0; i < size1; i++)
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
    delete []ut;
    delete []vx;
    delete []vy;
    delete []vt;

    delete []Ix;
    delete []Iy;
    delete []Iw;
    delete []Iwx;
    delete []Iwy;
    delete []Ixx;
    delete []Iyy;
    delete []Ixy;
    delete []Iwxx;
    delete []Iwyy;
    delete []Iwxy;

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
    delete []psi5;
    delete []psi6;
}


/**
  *
  * Function to normalize the images between 0 and 255
  *
**/
void image_normalization(
    const float *I,    //input images
    float       *In,   //normalized output images
    int          size  //size of the images
)
{
    //find the max and min values of the images
    const float max = *std::max_element(I, &I[size]);
    const float min = *std::min_element(I, &I[size]);
    const float den = max - min;

    if(den > 0)
	//normalize the images between 0 and 255
	for(int i = 0; i < size; i++)
	    In[i] = 255.0 * (I[i] - min) / den;

    else
	//copy the original data
	for(int i = 0; i < size; i++)
	    In[i] = I[i];
}


/**
  *
  *  Multiscale approach for computing the optical flow
  *
**/
void brox_optic_flow(
    const float *I,          //sequence of images
    float *u,                //x component of the optical flow
    float *v,                //y component of the optical flow
    const int    nxx,        //image width
    const int    nyy,        //image height
    const int    frames,     //number of frames
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
    if(frames <= 2)
    {
	std::cerr << "The method needs more than two frames" << std::endl;
	return;
    }

    const int size  = nxx * nyy * frames;

    std::vector<float *> Is(nscales);
    std::vector<float *> us(nscales);
    std::vector<float *> vs(nscales);

    std::vector<int> nx(nscales);
    std::vector<int> ny(nscales);

    Is[0] = new float[size];

    //normalize the input images between 0 and 255
    image_normalization(I, Is[0], size);

    //presmoothing the finest scale images
    for(int f = 0; f < frames; f++)
	gaussian(&Is[0][f*nxx*nyy], nxx, nyy, GAUSSIAN_SIGMA);

    us[0] = u;
    vs[0] = v;
    nx[0] = nxx;
    ny[0] = nyy;

    //create the scales
    for(int s = 1; s < nscales; s++)
    {
	zoom_size(nx[s-1], ny[s-1], nx[s], ny[s], nu);
	const int size  = nx[s] * ny[s] * frames;
	const int size1 = nx[s] * ny[s] * (frames-1);

	Is[s] = new float[size];
	us[s] = new float[size1];
	vs[s] = new float[size1];

	//compute the zoom from the previous scale
	for(int f = 0; f < frames; f++)
	    zoom_out(&Is[s-1][f*nx[s-1]*ny[s-1]], &Is[s][f*nx[s]*ny[s]], nx[s-1], ny[s-1], nu);
    }

    //initialization of the optical flow at the coarsest scale
    for(int i = 0; i < nx[nscales-1] * ny[nscales-1] * (frames-1); i++)
	us[nscales-1][i] = vs[nscales-1][i] = 0.0;

    //pyramidal approach for computing the optical flow
    for(int s = nscales-1; s >= 0; s--)
    {
	if(verbose) std::cout << "Scale: " << s << std::endl;
	
	//compute the optical flow for the current scale
	brox_optic_flow(
	    Is[s], us[s], vs[s], nx[s], ny[s], frames,
	    alpha, gamma, TOL, inner_iter, outer_iter, verbose
	);

	//if it is not the finer scale, then upsample the optical flow and adapt it conveniently
	if(s)
	{
	    for(int f = 0; f < frames-1; f++)
	    {
		zoom_in(&us[s][f*nx[s]*ny[s]], &us[s-1][f*nx[s-1]*ny[s-1]], nx[s], ny[s], nx[s-1], ny[s-1]);
		zoom_in(&vs[s][f*nx[s]*ny[s]], &vs[s-1][f*nx[s-1]*ny[s-1]], nx[s], ny[s], nx[s-1], ny[s-1]);
	    }

	    const int size = nx[s-1] * ny[s-1] * (frames-1);
	    
	    for(int i = 0; i < size; i++)
	    {
		us[s-1][i] *= 1.0 / nu;
		vs[s-1][i] *= 1.0 / nu;
	    }
	}
    }

    //delete allocated memory
    delete []Is[0];

    for(int i = 1; i < nscales; i++)
    {
	delete []Is[i];
	delete []us[i];
	delete []vs[i];
    }
}

#endif
