Robust Optical Flow
--------------------


*******
SUMMARY
*******

A program for optical flow estimation based on the Brox et al. method (spatial version):

Cite: Thomas Brox , Andrés Bruhn , Nils Papenberg , Joachim Weickert, "High Accuracy Optical Flow 
Estimation Based on a Theory for Warping", In Proc. 8th European Conference on Computer Vision, 
Springer LNCS 3024, T. Pajdla and J.Matas (Eds.), vol. 4, pp. 25-36, Prague, May 2004.

The program is part of an IPOL publication:
http://www.ipol.im/pub/algo/sms_optic_flow/

This program is written by 
Javier Sánchez Pérez <jsanchez@dis.ulpgc.es> CTIM, Universidad de Las Palmas de Gran Canaria 
Nelson Monzón López <nmonzon@ctim.es> CTIM, Universidad de Las Palmas de Gran Canaria
Agustín Salgado de la Nuez <asalgado@dis.ulpgc.es> CTIM, Universidad de Las Palmas de Gran Canaria

Version 1, released on June 20, 2012

This software is distributed under the terms of the BSD license (see file license.txt)



***********
COMPILATION
***********

Required environment: Any unix-like system with a standard compilation
environment (make and C and C++ compilers)

Required libraries: libpng, lipjpeg, libtiff

Compilation instructions: run "make" to produce an executable "main" 


*****
USAGE
*****

The program takes two input images, produce an optical flow as output, and
take some parameters.  The meaning of the parameters is thoroughly discussed on
the accompanying IPOL article.


Run:

	./main I1 I2 out_file

or

	./main I1 I2 out_file processors alpha gamma nscales zoom_factor TOL inner_iter outer_iter verbose

where:

	I1: first input image
	I2: second input image
	out_file: name of the output optical flow file
	processors: number of processors to run the method
	alpha: weight of the smoothing term
	gamma: weight of the gradient constancy term
	nscales: desired number of scales
	zoom_factor: downsampling factor 
	TOL: stopping criterion threshold for the numerical scheme
	inner_iter: number of inner iterations in the numerical scheme
	outer_iter: number of outer iterations in the numerical scheme
	verbose: 0 or 1, for quiet or verbose behaviour

Parameters can be ommited starting from the end of the list, and they will be
assigned reasonable default values.  Examples:

	./main I1.png I2.png flow.uv
	./main I1.png I2.png flow.uv 0 145 15 5 0.75 0.0001 1 38 1

If a parameter is given an invalid value it will take the default value.  If
the output filename is omitted, the flow will be saved on a file named
"flow.uv".
