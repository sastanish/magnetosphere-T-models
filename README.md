# Summary

This collection of files is a fortran wrapper of the magnetosphere field models of N. A. Tsyganenko. See his [website](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/) for the original fortran codes and documentation.  This project extends Tsyganenko's original codes with parallelized helper functions. These use the magnetosphere models to efficiently generate full-field data stored in standard netcdf format. Included with the original magnetosphere modelling functionallity, is a small program that can calculate the local reconnection rate as defined in [this](https://arxiv.org/abs/2502.01251) paper.

This directory contains:

 - *README.md* 
    - This file
 - *src/*
     - *main.f90*
       - The main program that feeds the models data and manages in/out and program control.
     - *inout.f90*
       - A simple module that implements input/output through the net-CDF library.
     - *TA16_RBF.f*
        - A modification of the original [TA16_RBF.f](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ta16/) model code. The compute block is wrapped in a module for interfacing with the main program. Redundant code (generating the RBF centers grid and loading parameters) have been placed in their own subroutines.
     - *TS05.for*
        - The original [TS04c.for](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ts05/) model code.  The only modification was to wrap the included subroutines inside a module interface.
     - *geopack.f*
        - As above, this is the original [geopack.f](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack) file. The only modification was to wrap the included subroutines inside a Module interface.
     - *recon.f90*
        - A short routine that reads in the netcdf field files and calculates the local reconnection rate.

# How to compile/install

This code is required to be compiled with [openMP](https://www.openmp.org/) and the fortran [NET-CDF library](https://docs.unidata.ucar.edu/netcdf-fortran/current/). Below are compilation commands for the Intel and GNU compilers.

## Intel
```
ifort inout.f90 geopack.f TA16_RBF.f TS05.f main.f90 -o TA_magneto.out -qopenmp `nf-config --fflags --flibs`
```
## Gnu Compilers
```
gfortran inout.f90 geopack.f TA16_RBF.f TS05.f main.f90 -o TA_magneto.out -fopenmp `nf-config --fflags --flibs`
```

Running the models require three files, 

  - *input_parameters.lst* -- The grid info for the models. 
  - *input_data.lst* -- The simulation input data prepared for the specific model being used (TS05 or TA16).
    - *TA16_RBF.par* -- The ta16 fitting parameters if the TA16 model is being used.

For an example of how to use this script, see the examples directory, examples/

To compile the reconnection code, the following work:
## Intel
```
ifort inout.f90 recon.f90 -o recon.out -qopenmp `nf-config --fflags --flibs`
```
## Gnu Compilers
```
gfortran inout.f90 recon.f90 -o recon.out -fopenmp `nf-config --fflags --flibs`
```


