#!/bin/bash
#======================================================
# Propogate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the standard partition (queue)
#SBATCH --partition=standard
#
# Specify project account
#SBATCH --account=mactaggart-aMHD
#
# Distribute processes in round-robin fashion
#SBATCH --distribution=cyclic
#
# Run the job on one node, all cores on the same node (full node)
#SBATCH --ntasks=1 --nodes=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=02:00:00
#
# Job name
#SBATCH --job-name=analysis
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

#./field_drop.out ../data/july_scans/s01/ 1 435
./plot_tail_slices.sh

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
