Horn-Schunck Optical Flow
-------------------------


*******
SUMMARY
*******

Two programs for optical flow estimation based on the Horn-Schunck functional.

These programs are part of an IPOL publication:
http://www.ipol.im/pub/algo/sm_horn_schunck/

Javier Sánchez Pérez <jsanchez@dis.ulpgc.es> CTIM, Universidad de Las Palmas de Gran Canaria
Enric Meinhardt-Llopis <enric.meinhardt@cmla.ens-cachan.fr> CMLA, ENS Cachan

Version 2, released on February 17, 2012

This software is distributed under the terms of the BSD license (see file license.txt)



***********
COMPILATION
***********

Required environment: Any unix-like system with a standard compilation
environment (make and C and C++ compilers)

Required libraries: libpng, libtiff

Compilation instructions: run "make" to produce two executables
"horn_schunck_classic" and "horn_schunck_pyramidal"


*****
USAGE
*****

All the programs take two input images, produce an optical flow as output, and
take some parameters.  The meaning of the parameters is thoroughly discussed on
the accompanying IPOL article.



1. Single scale Horn-Schunck.

Run:

	./horn_schunck_classic NITER ALPHA i0 i1 out

where

	NITER: number of iterations (e.g., 1000)
	ALPHA: weight of regularization term (e.g., 20)
	i0: first input image
	i1: second input image
	out: output optical flow (two-channel image)

All parameters are compulsory.  Example:

	./horn_schunck_classic 1000 30 a.png b.png flow.flo




2. Multiscale Horn-Schunck.

Run:

	./horn_schunck_pyramidal i0 i1 out

or

	./horn_schunck_pyramidal i0 i1 out NPROC ALPHA NSCALES ETA NWARPS EPSILON MAXITER verbose

where:

	i0: first input image
	i1: second input image
	out: output optical flow
	NPROC: number of OMP processors to run the method
	ALPHA: weight of the regularization term (e.g. 7)
	NSCALES: desired number of scales (e.g. 10)
	ETA: scale factor (e.g. 0.5)
	NWARPS: number of warps per scale (e.g. 10)
	EPSILON: tolerance of the stopping rule (e.g. 0.0001)
	MAXITER: maximum allowed number of iterations at each warp (e..g. 150)
	verbose: 0 or 1, for quiet or verbose behaviour

Parameters can be ommited starting from the end of the list, and they will be
assigned reasonable default values.  Examples:

	./horn_schunck_pyramidal a.png b.png flow.flo
	./horn_schunck_pyramidal a.png b.png flow.flo 0 100 5 0.5 50 0.0001 150 1

If a parameter is given a negative value it will take the default value.  If
the output filename is omitted, the flow will be saved on a file named
"flow.flo".
