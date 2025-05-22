#!/bin/bash
#======================================================
# Propogate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the standard partition (queue)
#SBATCH --partition=dev
#
# Specify project account
#SBATCH --account=mactaggart-aMHD1
#
# Distribute processes in round-robin fashion
#SBATCH --distribution=cyclic
#
# Run the job on one node, all cores on the same node (full node)
#SBATCH --ntasks=20 --nodes=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=00:15:00
#
# Job name
#SBATCH --job-name=calc_rad_rates
#
# Output file
#SBATCH --output=%j.out
#======================================================

module purge

#Example module load command. 
#Load any modules appropriate for your program's requirements

module load anaconda/python-3.10.9/2023.03

#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

conda activate magnetosphere

python track_x_point.py

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
