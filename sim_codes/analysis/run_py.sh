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
#SBATCH --ntasks=5 --nodes=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=20:00:00
#
# Job name
#SBATCH --job-name=plot_slices
#
# Output file
#SBATCH --output=%j-python.out
#======================================================

module purge

#Example module load command. 
#Load any modules appropriate for your program's requirements

module load anaconda/python-3.10.9/2023.03 intel/intel-2020.4 netcdf-fortran/intel-2020.4/4.5.4

#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

conda activate magnetosphere

python plot_tail_slices.py
#python calc_x-point_location.py
#python calc_total_nightside_field.py

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
