# KeplerMinMax
A novel closed-form solution to compute the surface altitude extrema for any two-body orbit is extended to oblate spheroids using Halley's method (a root finding technique with cubic convergence). Runtime performance has improved between three and five orders of magnitude while maintaining solution accuracy within several centimeters of numerically computed extrema when comparing this routine to conventional numerical minimization techniques such as fminbnd.  This method is simple to port to flight hardware it requires little computational overhead, few lines of code, and a small memory footprint.

This repository contains a c code implementation of keplerMinMax as outlined in the 2021 AMOS conference paper Novel Closed Form Solution for Orbit Segment Altitude Extrema Over Spherical and Oblate Central Bodies, it also has a MATLAB-C wrapper which allows for calling in a MATLAB environment. To compile a C MEX file in your MATLAB environment, execute the following command in your MATLAB command window:

    >> mex keplerMinMaxMex.c keplerMinMax.c
