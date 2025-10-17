# Compilation

This code is required to be compiled with the fortran NET-CDF library. Depending on your compiler of choice, the following command should take care of linking this code with it's dependencies as long as these are installed on your machine.

---
## Gnu

```
gfortran inout.f90 geopack.f TA16_RBF.f main.f90 -o main.out -fopenmp `nf-config --fflags --flibs`
```

## Intel
```
ifort inout.f90 geopack.f TA16_RBF.f main.f90 -o main.out -qopenmp `nf-config --fflags --flibs`
```
---

# Usage
The generated executable, ta16.out, requires three files;

  - *input_parameters.lst* -- The grid info for the ta16 model. 
  - *input_data.lst* -- The simulation input data prepared for the ta16 model
  - *TA16_RBF.par* -- The ta16 fitting parameters

For an example of how to use this script, see the examples directory, *$PROJECT_HOME/examples/*. A script to generate the correct ta16 input data from OMNI datasets is included in the *input_preparation* directory.
