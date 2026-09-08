# Gnu Compilers

This code is required to be compiled with open MP and the fortran NET-CDF library. The following command should take care of linking this code with it's dependencies as longs as these are installed on your machine.

```
gfortran inout.f90 geopack.f TA16_RBF.f main.f90 -o ta16.out -fopenmp `nf-config --fflags --flibs`
```

This generates a file called ta16.out. Running this requires three files, 

  - *input_parameters.lst* -- The grid info for the ta16 model. 
  - *input_data.lst* -- The simulation input data prepared for the ta16 model
  - *TA16_RBF.par* -- The ta16 fitting parameters

For an example of how to use this script, see the tests directory, *$PROJECT_HOME/sim_codes/tests*

The reconnection metrics are calculated in another file. To compile, run

```
gfortran inout.f90 recon.f90 -o recon.out -fopenmp `nf-config --fflags --flibs`

```
An example on how to run this program is in the tests directory.

## ARCHIE-W modules

On the Archie-West clusters, please ensure the netcdf-fortran libraries are loaded.

```
module load netcdf-fortran
```
