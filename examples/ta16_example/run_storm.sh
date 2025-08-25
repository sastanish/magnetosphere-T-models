#!/bin/bash
#======================================================
# Propogate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the standard partition (queue)
#SBATCH --partition=dev
#
# Specify project account
#SBATCH --account=mactaggart-aMHD
#
# Distribute processes in round-robin fashion
#SBATCH --distribution=cyclic
#
# Run the job on one node, all cores on the same node (full node)
#SBATCH --ntasks=20 --nodes=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=00:20:00
#
# Job name
#SBATCH --job-name=storm_test
#
# Output file
#SBATCH --output=%j_cmd.out
#======================================================

module purge

#Example module load command. 
#Load any modules appropriate for your program's requirements

module load netcdf-fortran/intel-2020.4 intel/intel-2020.4

#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

export OMP_NUM_THREADS=20
./ta16.out 3 10

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
