# Summary

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
     - *TS04c.for*
        - The original [TS04c.for](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ts05/) model code.  The only modification was  to wrap the included subroutines inside a Module interface.
     - *geopack.f*
        - As above, this is the original [geopack.f](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack) file. The only modification was to wrap the included subroutines inside a Module interface.
    - *recon.f90*
        - This file contains fortran routines for calculating the reconnection rate as written in [this paper](https://arxiv.org/abs/2502.01251).

# How to compile/install

This code is required to be compiled with open MP and the fortran NET-CDF library. The following command should take care of linking this code with it's dependencies as long as these are installed on your machine.

## Intel
```
ifort inout.f90 geopack.f TA16_RBF.f main.f90 -o ta16.out -qopenmp `nf-config --fflags --flibs`
```
## Gnu Compilers
```
gfortran inout.f90 geopack.f TA16_RBF.f main.f90 -o ta16.out -fopenmp `nf-config --fflags --flibs`
```

This generates a file called ta16.out. Running this requires three files, 

  - *input_parameters.lst* -- The grid info for the ta16 model. 
  - *input_data.lst* -- The simulation input data prepared for the ta16 model
  - *TA16_RBF.par* -- The ta16 fitting parameters

For an example of how to use this script, see the examples directory, examples/

The reconnection metrics are calculated in another file. To compile, run

```
ifort inout.f90 recon.f90 -o recon.out -qopenmp `nf-config --fflags --flibs`

```
An example on how to run this program is in the examples directory.
