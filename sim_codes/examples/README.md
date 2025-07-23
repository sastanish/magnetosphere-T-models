This tutorial will walk you through how to use the TA16 model and attached reconnection metrics program. First, compile both the *recon.f90* and *main.f90* programs in the *src/* directory and copy their binaries to this directory. These programs will be referred to as *ta16.out* and *recon.out*. To run these programs, 

  - 1) Write an *input_parameters.txt* file which contains the grid info. See the example here.

  - 2) Supply an *input_data.lst* as generated for the ta16 model

  - 3) Ensure that the *TA16_RBF.par* parameter file is in the same directory you will run the scripts from

  - 4) Set the number of cores being used by exporting the variable **OMP_NUM_THREADS=**. 

  - 5) Run the *ta16.out* program. This will generate a series of netcdf files named *output_1.nc*, *output_2.nc*, etc. Each of these contains the 3 components of the magnetic field. You must specify the start and end line numbers for the TA16 model from the command line. Input two integers for the line to start and end the computation on in the *input_data.lst* file. For example,
  ```
  ./ta16.out 1 100
  ```
  This will execute the ta16 model on lines 1-100 and output files *output_1.nc* - *output_100.nc*

  - 6) Run the *recon.out* program specifying which output files we with to calculate. For example, 
    ```
    ./recon.out 1 100
    ```
    will calculate the reconnection metrics over the output files numbered *output_1.nc* up to *output_100.nc*. These metrics will be saved in files named *rate_1.nc*, etc.

An example slurm script for the Archie-West machine is included as *run_storm.sh*.
