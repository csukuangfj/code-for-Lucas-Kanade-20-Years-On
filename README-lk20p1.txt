Lucas-Kanade 20 Years On: A Unifying Framework

Part 1: The Quantity Approximated, the Warp Update Rule and the Gradient
Descent Approximation

Simon Baker and Iain Matthews

The Robotics Institute
Carnegie Mellon University
Pittsburgh, PA 15213
USA

Contact: simonb@cs.cmu.edu, iainm@cs.cmu.edu

---------------------------------------------------------------------------

URL: http://www.ri.cmu.edu/projects/project_515.html

Bibtex reference:

@article{Baker_2004_4293,
	author = "Simon Baker and Iain Matthews",
	title = "Lucas-Kanade 20 Years On: A Unifying Framework Part 1: The Quantity Approximated, the Warp Update Rule, and the Gradient Descent Approximation",
	journal = "International Journal of Computer Vision",
	year = "2004"
}


---------------------------------------------------------------------------
                             MATLAB SOURCE CODE
---------------------------------------------------------------------------

The source code is split into two directories, both need to be added to your
Matlab search path:

lk20-p1           The algorithms described in the paper and test scripts
lk20-common       Common code for all algorithms (and for later papers)


lk20-p1/
--------

The Gauss-Newton affine warp algorithms in Section 3 of the paper are:

affine_fa.m       Forwards-Additive
affine_fc.m       Forwards-Compositional
affine_ia.m       Inverse-Additive
affine_ic.m       Inverse-Compositional

The Gauss-Newton projecive warp algorithms in Section 3 of the paper are:

homo_fa.m         Forwards-Additive
homo_fc.m         Forwards-Compositional
homo_ic.m         Inverse-Compositional

Note: the inverse-additive algorithm cannot be applied with a projective warp.


The additional algorithms from Section 4 of the paper are:

affine_ic_nt.m    Inverse-Compositional Newton (4.2)
affine_ic_sd.m    Inverse-Compositional Steepest Descent (4.3)
affine_ic_d.m     Inverse-Compositional with diagonal Hessian (4.4)
affine_ic.nt_d.m  Inverse-Compositional Newton with diagonal Hessian (4.4)
affine_ic_lm.m    Inverse-Compositional Levenburg-Marquardt (4.5)


To run the perturbation experiments described in the paper examine:

run_affine.m       Run an affine warp experiment
run_homo.m         Run a projective warp experiment

These files should be edited for the experimental parameters you would like to
test. By default they run everything for a few tests. To replicate the
experiments in the paper run each algorithm for 100 convergence tests and 1000
frequency of convergence tests. This will take some time!


lk20p1/data/
------------

The data directory contains precomputed random perturbation data (for
experimental consistency) and the images and template definitions.


lk20p1/figures/
---------------

The figures directory contains matlab scripts for plotting your results and
creating graphs similar to those shown in the paper.
