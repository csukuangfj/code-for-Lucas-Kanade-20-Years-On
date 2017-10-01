IPOL 2D CLG Optical Flow Demo Program

This program implements the combined local-global (CLG) algorithm for optical 
flow computation in 2D greyscale images, introduced by Bruhn et al. (2002). 
The code is supplied in addition to the IPOL documentation material, available 
online at http://www.ipol.im.

- Program written by Jorge Jara <jjara@dcc.uchile.cl>, José Delpiano 
<jdelpian@uandes.cl> and Mauricio Cerda <mauriciocerda@med.uchile.cl>. 
Distributed under the terms of the BSD license. See the file license.txt 
for details.

- The files iio.h and iio.c are Copyright (c) 2012, Enric Meinhardt Llopis. 
Distributed under the terms of the BSD license. See the file license.txt 
for details.

- The files bicubic_interpolation.c, bicubic_interpolation.h, mask.c, mask.h, 
zoom.c, and zoom.h are Copyright (c) 2012, Javier Sánchez.
Distributed under the terms of the BSD license. See the file license.txt 
for details.


== Contents ==

bin             executable folder (empty, intended to be used with CMake)
img             sample images folder
source          source code folder

CMakeLists.txt  CMake configuration file (there is also a configuration file in 
                the source folder)
iio_license.txt license file for the utility library iio
license.txt     license file for the CLG optical flow program


== Usage ==

From command line:

clg_of <image1> <image2> <alpha> <rho> <sigma> <iterations> 
    [<w> <n_scales> <zoom_factor> <coupledMode> <verbose>] <flow_output_image>

Note that if the program was compiled using the provided CMake script, the 
executable will be located in the "bin" folder.

Parameters:

<image1>, <image2> specify the input time frames (images) for the optical flow 
                   computation.

<alpha>            specifies the value of the global regularization coefficient.
                   Use a value greater than 0.

<rho>              specifies the value of the local spatio-temporal derivatives 
                   regularization coefficient in the form of a gaussian 
                   smoothing with standard deviation rho (greater than 0). Use 
                   0 to avoid this smoothing.

<sigma>            specifies the standard deviation value (greater than 0) for a 
                   gaussian smoothing applied to the input images. Use 0 to 
                   avoid smoothing.

<iterations>       max. number of iterations for the optical flow computation.

<w>                SOR relaxation factor, between 0 and 2. Default: 1.8.

<n_scales>         mumber of scales for the pyramidal approach. Default: 1.
                   The number is verified to have at most a size 30 pixels
                   image at the highest pyramid level.

<zoom_factor>      scale factor, between 0 and 1. Default: 0.5.

<coupledMode>      iteration type. 1->PCGS, 0->SOR. Default: 1.

<verbose>          shows (1) or hide (0) messages. Default: 1.
                   Also, the output is saved in two text files (for scripting)
                   and can be readed in Matlab with:
                       vx=load(file1, '-ascii')';
                       vy=load(file2, '-ascii')';

<flow_output_image> specifies a TIFF output image containing the output flow 
                    magnitudes.

If needed, first convert RGB images to grayscale with the following commands:

test_rgb2gray <image1> <image1_gray>
test_rgb2gray <image2> <image2_gray>

When using the CMake script the excutable will be located in the "bin" folder.


== Use examples ==

Basic example:
./bin/test_clgof img/baboon_rotation/a_bw.png img/baboon_rotation/b_bw.png \
200 5.0 1.0 200 flow.tiff

In this example alpha=200, rho=5.0, sigma=1.0, iteartions=200, and output saved 
in flow.tiff

Full example:
./bin/test_clgof img/baboon_rotation/a_bw.png img/baboon_rotation/b_bw.png \
200 5.0 1.0 300 1.0 2 .65 1 1 flow.tiff

In this example alpha=200, rho=5.0, sigma=1.0, iteartions=200, and output saved 
in flow.tiff. Also w = 1.0 (Gauss-Seidel iteration), with 2 scales, 
downsampling factor of .65, using coupled Gaussian iteration (w is ignored), 
and verbose mode activated.


== Compiling ==

Instructions are included below for compiling on Linux sytems with GCC, on 
Windows with MSVC.


== Compiling (CMake) ==

The  program can be compiled using CMake and make. CMake can be obtained at

http://www.cmake.org

The compilation can be performed by entering the following commands:

cd build 
cmake ..
make

Other  compilers and libraries required are described below.

== Compiling (Linux) ==

To compile this software under Linux, first install the development files for
libjpeg, libpng, and libtiff.  On Ubuntu and other Debian-based systems, enter
the following into a terminal:
    sudo apt-get install build-essential libjpeg8-dev libpng-dev libtiff-dev
On Redhat, Fedora, and CentOS, use
    sudo yum install make gcc libjpeg-turbo-devel libpng-devel libtiff-devel

--- Building with PNG, JPEG, and/or TIFF support ---

To use the program with PNG, JPEG, and/or TIFF images, the 
following libraries are needed.

    For PNG:    libpng and zlib
    For JPEG:   libjpeg 
    For TIFF:   libtiff

These libraries can be obtained at 
    
    http://www.libpng.org/pub/png/libpng.html
    http://www.zlib.net/
    http://www.ijg.org/
    http://www.remotesensing.org/libtiff/


== Compiling (Windows with MSVC) ==

The express version of the Microsoft Visual C++ (MSVC) compiler can be 
obtained for free at

    http://www.microsoft.com/visualstudio/en-us/products/2010-editions/express

To use the program with PNG, JPEG, and/or TIFF images, the 
following libraries are needed.

    For PNG :  libpng and zlib
    For JPEG:  libjpeg 
    For TIFF:  libtiff

You will have to specify the location of the corresponding libraries when 
compiling.

== Acknowledgements ==

This material is based upon work supported by Comisión Nacional de 
Investigación Científica y Tecnológica (CONICYT), Chile, through FONDECYT 
projects No. 1090246, 1120579, 3140447.
Jorge Jara and José Delpiano funded by CONICYT PhD scholarships.
Any opinions, findings, and conclusions or recommendations expressed in this 
material are those of the author(s) and do not necessarily reflect the views 
of the CONICYT.
