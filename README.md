# KeplerMinMax
A novel closed-form solution to compute the surface altitude extrema for any two-body orbit is extended to oblate spheroids using Halley's method (a root finding technique with cubic convergence).  I've noticed a runtime performance improvement between three and five orders of magnitude while maintaining solution accuracy within several centimeters of numerically computed extrema.  This method is simple to port to flight hardware it requires little computational overhead, few lines of code, and a small memory footprint.

This repository contains a c code implementation of keplerMinMax as outlined in the 2021 AMOS conference paper Novel Closed Form Solution for Orbit Segment Altitude Extrema Over Spherical and Oblate Central Bodies.
