// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#ifndef MASK_H
#define MASK_H

/**
 *
 * Function to apply a 3x3 mask to an image
 *
 */
void mask3x3(
    const float *input, //input image
    float *output,      //output image
    const int nx,       //image width
    const int ny,       //image height
    const float *mask   //mask to be applied
)
{
    //apply the mask to the center body of the image
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    double sum = 0;
	    for(int l = 0; l < 3; l++)
	    {
		for(int m = 0; m < 3; m++)
		{
		    int p = (i + l -1) * nx + j + m -1;
		    sum += input[p] * mask[l * 3 + m];
		}
	    }
	    int k = i * nx + j;
	    output[k] = sum;
	}
    }

    //apply the mask to the first and last rows
    for(int j = 1; j < nx-1; j++)
    {
	double sum = 0;
	sum += input[j-1] * (mask[0] + mask[3]);
	sum += input[ j ] * (mask[1] + mask[4]);
	sum += input[j+1] * (mask[2] + mask[5]);

	sum += input[nx + j-1] * mask[6];
	sum += input[nx +  j ] * mask[7];
	sum += input[nx + j+1] * mask[8];

	output[j] = sum;

	sum = 0;
	sum += input[(ny-2)*nx+j-1] * mask[0];
	sum += input[(ny-2)*nx+j  ] * mask[1];
	sum += input[(ny-2)*nx+j+1] * mask[2];

	sum += input[(ny-1)*nx+j-1] * (mask[6] + mask[3]);
	sum += input[(ny-1)*nx+j  ] * (mask[7] + mask[4]);
	sum += input[(ny-1)*nx+j+1] * (mask[8] + mask[5]);

	output[(ny-1)*nx + j] = sum;
    }

    //apply the mask to the first and last columns
    for(int i = 1; i < ny-1; i++)
    {
	double sum = 0;
	sum += input[(i - 1)*nx]   * (mask[0] + mask[1]);
	sum += input[(i - 1)*nx+1] * mask[2];

	sum += input[i * nx]   * (mask[3] + mask[4]);
	sum += input[i * nx+1] * mask[5];

	sum += input[(i + 1)*nx]   * (mask[6] + mask[7]);
	sum += input[(i + 1)*nx+1] * mask[8];

	output[i*nx] = sum;

	sum = 0;
	sum += input[i * nx-2] * mask[0];
	sum += input[i * nx-1] * (mask[1] + mask[2]);

	sum += input[(i + 1)*nx-2] * mask[3];
	sum += input[(i + 1)*nx-1] * (mask[4] + mask[5]);

	sum += input[(i + 2)*nx-2] * mask[6];
	sum += input[(i + 2)*nx-1] * (mask[7] + mask[8]);

	output[i*nx + nx -1] = sum;
    }

    //apply the mask to the four corners
    output[0] = input[0]    * (mask[0] + mask[1] + mask[3] + mask[4]) +
	input[1]    * (mask[2] + mask[5]) +
	input[nx]   * (mask[6] + mask[7]) +
	input[nx+1] * mask[8];

    output[nx-1] =
	input[nx-2]   * (mask[0] + mask[3]) +
	input[nx-1]   * (mask[1] + mask[2] + mask[4] + mask[5]) +
	input[2*nx-2] * mask[6] +
	input[2*nx-1] * (mask[7] + mask[8]);

    output[(ny-1)*nx] =
	input[(ny-2)*nx]   * (mask[0] + mask[1]) +
	input[(ny-2)*nx+1] *  mask[2] +
	input[(ny-1)*nx]   * (mask[3] + mask[4] + mask[6] + mask[7]) +
	input[(ny-1)*nx+1] * (mask[5] + mask[8]);

    output[ny*nx-1] =
	input[(ny-1)*nx-2] * mask[0] +
	input[(ny-1)*nx-1] * (mask[1] + mask[2]) +
	input[ny*nx-2] * (mask[3] + mask[6]) +
	input[ny*nx-1] * (mask[4] + mask[5] + mask[7] + mask[8]);
}

/**
 *
 * Compute the second order X derivative
 *
 */
void Dxx(
    const float *I, //input image
    float *Ixx,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
)
{
    //mask of second derivative
    float M[]  = {0., 0., 0.,
		  1.,-2., 1.,
		  0., 0., 0.};

    //computing the second derivative
    mask3x3(I, Ixx, nx, ny, M);
}


/**
 *
 * Compute the second order Y derivative
 *
 */
void Dyy(
    const float *I, //input image
    float *Iyy,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
)
{
    //mask of second derivative
    float M[]  = {0., 1., 0.,
		  0.,-2., 0.,
		  0., 1., 0.};

    //computing the second derivative
    mask3x3(I, Iyy, nx, ny, M);
}


/**
 *
 * Compute the second order XY derivative
 *
 */
void Dxy(
    const float *I, //input image
    float *Ixy,     //oputput derivative
    const int nx,   //image width
    const int ny    //image height
)
{
    //mask of second derivative
    float M[]  = {1./4., 0.,-1./4.,
		  0.,    0., 0.,
		  -1./4., 0., 1./4.};

    //computing the second derivative
    mask3x3(I, Ixy, nx, ny, M);
}


/**
 *
 * Compute the gradient with central differences
 *
 */
void gradient(
    const float *input, //input image
    float *dx,          //computed x derivative
    float *dy,          //computed y derivative
    const int nx,       //image width
    const int ny        //image height
)
{
    //gradient in the center body of the image
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    const int k = i * nx + j;
	    dx[k] = 0.5*(input[k+1] - input[k-1]);
	    dy[k] = 0.5*(input[k+nx] - input[k-nx]);
	}
    }

    //gradient in the first and last rows
    for(int j = 1; j < nx-1; j++)
    {
	dx[j] = 0.5*(input[j+1] - input[j-1]);
	dy[j] = 0.5*(input[j+nx] - input[j]);

	const int k = (ny - 1) * nx + j;

	dx[k] = 0.5*(input[k+1] - input[k-1]);
	dy[k] = 0.5*(input[k] - input[k-nx]);
    }

    //gradient in the first and last columns
    for(int i = 1; i < ny-1; i++)
    {
	const int p = i * nx;
	dx[p] = 0.5*(input[p+1] - input[p]);
	dy[p] = 0.5*(input[p+nx] - input[p-nx]);

	const int k = (i+1) * nx - 1;

	dx[k] = 0.5*(input[k] - input[k-1]);
	dy[k] = 0.5*(input[k+nx] - input[k-nx]);
    }

    //calculate the gradient in the corners
    dx[0] = 0.5*(input[1] - input[0]);
    dy[0] = 0.5*(input[nx] - input[0]);

    dx[nx-1] = 0.5*(input[nx-1] - input[nx-2]);
    dy[nx-1] = 0.5*(input[2*nx-1] - input[nx-1]);

    dx[(ny-1)*nx] = 0.5*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
    dy[(ny-1)*nx] = 0.5*(input[(ny-1)*nx] - input[(ny-2)*nx]);

    dx[ny*nx-1] = 0.5*(input[ny*nx-1] - input[ny*nx-1-1]);
    dy[ny*nx-1] = 0.5*(input[ny*nx-1] - input[(ny-1)*nx-1]);
}


/**
 *
 * Compute the coefficients of the divergence term
 *
 */
void psi_divergence(
    const float *psi, //robust functional
    float *psi1,      //coefficients of divergence
    float *psi2,      //coefficients of divergence
    float *psi3,      //coefficients of divergence
    float *psi4,      //coefficients of divergence
    const int nx,     //image width
    const int ny      //image height
)
{
    //calculate coefficients in the center body of the image
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    const int k = i * nx + j;

	    psi1[k] = 0.5 * (psi[k +nx] + psi[k]);
	    psi2[k] = 0.5 * (psi[k -nx] + psi[k]);
	    psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
	    psi4[k] = 0.5 * (psi[k - 1] + psi[k]);
	}
    }

    //calculate coefficients in the first and last rows
    for(int j = 1; j < nx-1; j++)
    {
	psi1[j] = 0.5 * (psi[j +nx] + psi[j]);
	psi2[j] = 0;
	psi3[j] = 0.5 * (psi[j + 1] + psi[j]);
	psi4[j] = 0.5 * (psi[j - 1] + psi[j]);

	const int k  = (ny-1)*nx + j;
	
	psi1[k] = 0;
	psi2[k] = 0.5 * (psi[k -nx] + psi[k]);
	psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
	psi4[k] = 0.5 * (psi[k - 1] + psi[k]);
    }

    //calculate coefficients in the first and last columns
    for(int i = 1; i < ny-1; i++)
    {
	const int k  = i*nx;
	
	psi1[k] = 0.5 * (psi[k +nx] + psi[k]);
	psi2[k] = 0.5 * (psi[k -nx] + psi[k]);
	psi3[k] = 0.5 * (psi[k + 1] + psi[k]);
	psi4[k] = 0;

	const int   j  = (i+1) * nx - 1;

	psi1[j] = 0.5 * (psi[j +nx] + psi[j]);
	psi2[j] = 0.5 * (psi[j -nx] + psi[j]);
	psi3[j] = 0;
	psi4[j] = 0.5 * (psi[j - 1] + psi[j]);
    }

    //up-left corner (0,0)
    psi1[0] = 0.5 * (psi[nx] + psi[0]);
    psi3[0] = 0.5 * (psi[1] + psi[0]);
    psi2[0] = psi4[0] = 0;

    //up-right corner (nx,0)
    psi1[nx-1] = 0.5 * (psi[nx-1 + nx] + psi[nx-1]);
    psi4[nx-1] = 0.5 * (psi[nx - 2] + psi[nx-1]);
    psi2[nx-1] = psi3[nx-1] = 0;

    //bottom-left corner (0,ny)
    psi2[(ny-1)*nx] = 0.5 * (psi[(ny-2) * nx] + psi[(ny-1) * nx]);
    psi3[(ny-1)*nx] = 0.5 * (psi[(ny-1)*nx + 1] + psi[(ny-1)*nx]);
    psi1[(ny-1)*nx] = psi4[(ny-1)*nx] = 0;

	    
    //bottom-right corner (nx,ny)
    psi2[ny*nx-1] = 0.5 * (psi[ny * nx - 1 - nx] + psi[ny * nx -1]);
    psi4[ny*nx-1] = 0.5 * (psi[ny * nx - 2] + psi[ny * nx -1]);
    psi1[ny*nx-1] = psi3[ny*nx-1] = 0;
}



/**
 *
 * Compute the divergence of the optical flow
 *
 */
void divergence_u(
    const float *u,    //x component of optical flow
    const float *v,    //y component of optical flow
    const float *psi1, //coefficients of divergence
    const float *psi2, //coefficients of divergence
    const float *psi3, //coefficients of divergence
    const float *psi4, //coefficients of divergence
    float *div_u,      //computed divergence for u
    float *div_v,      //computed divergence for v
    const int nx,      //image width 
    const int ny       //image height
)
{
    //calculate the divergence in the center body of the image
    for(int i = 1; i < ny-1; i++)
    {
	for(int j = 1; j < nx-1; j++)
	{
	    const int k = i * nx + j;

	    div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) + 
			psi3[k] * (u[k + 1]  - u[k]) + psi4[k] * (u[k - 1]  - u[k]);
	    div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) + 
			psi3[k] * (v[k + 1]  - v[k]) + psi4[k] * (v[k - 1]  - v[k]);
	}
    }

    //calculate the divergence in the first and last rows
    for(int j = 1; j < nx-1; j++)
    {
	div_u[j] = psi1[j] * (u[j + nx] - u[j]) + psi3[j] * (u[j + 1]  - u[j]) + psi4[j] * (u[j - 1] - u[j]);
	div_v[j] = psi1[j] * (v[j + nx] - v[j]) + psi3[j] * (v[j + 1]  - v[j]) + psi4[j] * (v[j - 1] - v[j]);

	const int   k  = (ny-1)*nx + j;

	div_u[k] = psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]) + psi4[k] * (u[k - 1] - u[k]);
	div_v[k] = psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]) + psi4[k] * (v[k - 1] - v[k]);	
    }

    //calculate the divergence in the first and last columns
    for(int i = 1; i < ny-1; i++)
    {
	const int   k  = i*nx;

	div_u[k] = psi1[k] * (u[k + nx] - u[k]) + psi2[k] * (u[k - nx] - u[k]) + psi3[k] * (u[k + 1] - u[k]);
	div_v[k] = psi1[k] * (v[k + nx] - v[k]) + psi2[k] * (v[k - nx] - v[k]) + psi3[k] * (v[k + 1] - v[k]);

	const int   j  = (i+1) * nx - 1;

	div_u[j] = psi1[j] * (u[j + nx] - u[j]) + psi2[j] * (u[j - nx] - u[j]) + psi4[j] * (u[j - 1] - u[j]);
	div_v[j] = psi1[j] * (v[j + nx] - v[j]) + psi2[j] * (v[j - nx] - v[j]) + psi4[j] * (v[j - 1] - v[j]);
    }

    //up-left corner (0,0)
    div_u[0] = psi1[0] * (u[nx] - u[0]) + psi3[0] * (u[1] - u[0]);
    div_v[0] = psi1[0] * (v[nx] - v[0]) + psi3[0] * (v[1] - v[0]);

    //up-right corner (nx,0)
    div_u[nx-1] = psi1[nx-1] * (u[nx - 1 + nx] - u[nx-1]) + psi4[nx-1] * (u[nx - 2]  - u[nx - 1]);
    div_v[nx-1] = psi1[nx-1] * (v[nx - 1 + nx] - v[nx-1]) + psi4[nx-1] * (v[nx - 2]  - v[nx - 1]);

    //bottom-left corner (0,ny)
    div_u[(ny-1)*nx] = psi2[(ny-1)*nx] * (u[(ny-2) * nx]     - u[(ny-1) * nx]) + 
			psi3[(ny-1)*nx] * (u[(ny-1) * nx + 1] - u[(ny-1) * nx]);
    div_v[(ny-1)*nx] = psi2[(ny-1)*nx] * (v[(ny-2) * nx]     - v[(ny-1) * nx]) + 
			psi3[(ny-1)*nx] * (v[(ny-1) * nx + 1] - v[(ny-1) * nx]);
	    
    //bottom-right corner (nx,ny)
    div_u[ny*nx-1] = psi2[ny*nx-1] * (u[ny * nx - 1 - nx] - u[ny * nx - 1]) + 
		      psi4[ny*nx-1] * (u[ny * nx - 2]      - u[ny * nx - 1]);
    div_v[ny*nx-1] = psi2[ny*nx-1] * (v[ny * nx - 1 - nx] - v[ny * nx - 1]) + 
		      psi4[ny*nx-1] * (v[ny * nx - 2]      - v[ny * nx - 1]);
}

#endif
