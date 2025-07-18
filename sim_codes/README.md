## Summary

This collection of files is a fortran wrapper of the TS06 solar field model and dipole codes by N. A. Tsyganenko. See his [website](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/) for the original fortran codes and documentation.  This directory contains:

 - *README.md* 
    - This file
 - *src/*
     - *main.f90*
       - The main program that feeds the TA16 model it's data and manages in/out and program control.
     - *inout.f90*
       - A simple module that implements input/output through the net-CDF library.
     - *TA16_RBF.f*
        - A modification of the original [TA16_RBF.f](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ta16/) model code. The compute block is wrapped in a module for interfacing with the main program. Redundant code (generating the RBF centers grid and loading parameters) have been placed in their own subroutines.
     - *geopack.f*
        - As above, this is the original [geopack.f](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack) file. The only modification was to wrap the included subroutines inside a Module interface.
    - *recon.f90*
        - This file contains fortran routines for calculating the reconnection rate as written in [this paper](https://arxiv.org/abs/2502.01251).
