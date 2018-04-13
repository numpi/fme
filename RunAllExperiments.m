% This script run all the experiments for the paper "Fast solvers for
% 2D fractional diffusion equations using rank structured matrices".
%
%

% Example 1: Testing the 1D case. This will produce files e1.dat and e2.dat
Experiment1D;

% Example 2: Testing the FD for the 2D problem, constant coefficient case.
% This will produce fd-time*.dat files
FD_Example;

% Example 3: The same as 2, in the variable coefficient case. This will
% produce fd-times-vc*.dat files.
FD_Example_vc;

% Example 4: Finite elements, essentially the same configuration of
% Example 2, just with a different spatial discretization. This will
% produce fe-times*.dat files.
FE_Example;
