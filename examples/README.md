## Summary

This collection of files is an example of using the python wrapper of the TS05/TA16 solar field models and dipole codes by N. A. Tsyganenko. See his [website](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/) for the original fortran codes and documentation.  This directory contains:

 - *README.md* 
    - This file
 - *example_may_storm_TS05/*
     - *may_2024_storm_with_TS05_vars.dat*
        - May 2024 storm for example data.
     - *storm.py*
        - An example calculation using the Tasyganenko Fortran wrapper for TS05.  
 - *example_may_storm_TA16/*
     - *convert_TS05_data_to_TA16_data.py*
        - Generates the input parameters for the TA16 model from a TS05 input file.
    - *may_2024_storm_with_TA16_vars.dat*
        - May 2024 storm example data
    - *TA16_RBF.par*
        - Model parameters for TA16 model. Must be present in working directory from where the model is run.
    - *storm.py*
        - An example calculation using the Tasyganenko Fortran wrapper for TA16.  
 - *fortran_example_storm.py*
    - An example calculation using the original Tasyganenko codes in modern fortran. Slower than then wrapped python (at the moment).

## Requirements

This example requires the following python packages:

 - xarray
 - multiprocess
 - pandas
 - netCDF4
 - h5netcdf - for compression when writing the magnetic field to disk.
 - The Tsyganenko_wrapper package from the modules_dir
