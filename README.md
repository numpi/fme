# Solving fractional PDEs by rank-structured matrices

This repository contains the MATLAB code required to reproduce the numerical experiments contained in the paper "*Fast solvers for 2D fractional diffusion equations using rank structured matrices*", S. Massei, M. Mazza, L. Robol, 2018. 

## Code structure
To reproduce the tests, you will need to download and add to your MATLAB path the [`hm-toolbox`](https://github.com/numpi/hm-toolbox) and the [`rktoolbox`](http://guettel.com/rktoolbox/) packages. Installation of those is as simple as unpacking them somewhere, and adding them to your MATLAB path by:

    >> addpath /path/to/hm-toolbox;
    >> addpath /path/to/rktoolbox; addpath /path/to/rktoolbox/utils;

Once those are loaded, you can run the various tests by running the scripts listed below: 

 1. `Experiment1D.m` replicates the 1D tests that compare the HODLR solver with the preconditioned GMRES. Two files `e1.dat` and `e2.dat` will be generated containing the data included in the table in the paper. 
 2. `FD_Example.m` and `FD_Example_vc.m` contain the examples of the 2D solver for the finite difference formulation, discretized using implicit Euler in time. The files implement the constant and variable coefficients case, respectively.
 3. `FE_Example.m` contains the finite element case, and is completely analogous to `FD_Example.m`. 

You might want to inspect the file <code>RunAllExperiments.m</code> for further information. 

## Contact
Did you find bugs, or have any kind of feedback? Let us know! 

[stefano.massei@epfl.ch](mailto:stefano.massei@epfl.ch), [mariarosa.mazza@ipp.mpg.de](mailto:mariarosa.mazza@ipp.mpg.de), [leonardo.robol@isti.cnr.it](mailto:leonardo.robol@isti.cnr.it)


