## Summary

This collection of files is a python wrapper of the TS05 solar field model and dipole codes by N. A. Tsyganenko. See his [website](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/) for the original fortran codes and documentation.  This directory contains:

 - *README.md* 
    - This file
 - *src/*
     - *TsyganenkoWrapper.f90*
        - A fortran wrapper written in F90 that serves as an interface between Tsyganenko's original F77 codes and python. This file has been written to accept an array of spacial values for which to calculate the magnetic field from. This allows us to avoid looping in python (slow!) without modifying Tsyganenko's code an any way.
     - *TsyganenkoWrapper.pyf*
        - This file serves as an interpreter between python and fortran, providing information about the KIND and INTENT of each variable in the *Tsyganenko_wrapper*.  
     - *TS04c.for*
        - The original [TS04c.for](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ts05/) model code.  The only modification was  to wrap the included subroutines inside a Module interface.
     - *geopack.f*
        - As above, this is the original [geopack.f](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack) file. The only modification was to wrap the included subroutines inside a Module interface.

## How to compile/install

To use these subroutines in python, all the fortran codes need compiled via the numpy wrapper, [f2py](https://numpy.org/doc/stable/f2py/f2py.getting-started.html).  You will need an active python environment with the following packages:

 - numpy
 - meson
 - ninja

You will also need a fortran compiler.  Ensure that you python environment is active and run,
```
python -m numpy.f2py -c TS04c.for TA16_RBF.f geopack.f reconnectionMetrics.f90 magnetosphereModels.f90 TsyganenkoWrapper.pyf TsyganenkoWrapper.f90
```
within the src directory.  This should generate a module file entitled *TsyganenkoWrapper.cpython-313-x86_64-linux-gnu.so* or something similar. This module file can now be imported into python. For how to use, see the example storm directory.

**NOTE:** Do not move this file around. If you wish to call this package in another directory, use a system link to point to this file in place.  Ie, *ln -s PATH_TO_MODULE CURRENT_DIR/.*
