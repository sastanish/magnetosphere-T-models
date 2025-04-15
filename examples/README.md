## Summary

This collection of files is an example of using the python wrapper of the TS05 solar field model and dipole codes by N. A. Tsyganenko. See his [website](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/) for the original fortran codes and documentation.  This directory contains:

 - *README.md* 
    - This file
 - *may_2024_storm_with_TS05_vars.dat*
    - May 2024 storm for example data.
 - *python_example_storm.py*
    - An example calculation using the Tasyganenko Fortran wrapper.  
 - *fortran_example_storm.py*
    - An example calculation using the original Tasyganenko codes in modern fortran.
## Requirements

This example requires the following python packages:

 - xarray
 - multiprocess
 - pandas
 - netCDF4
 - h5netcdf - for compression when writing the magnetic field to disk.
 - The Tsyganenko_wrapper package from the modules_dir
