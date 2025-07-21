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
#SBATCH --ntasks=1 --nodes=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=00:30:00
#
# Job name
#SBATCH --job-name=plot_slices
#
# Output file
#SBATCH --output=%j_cmd.out
#======================================================

module purge

#Example module load command. 
#Load any modules appropriate for your program's requirements

module load netcdf-fortran/intel-2020.4 intel/intel-2020.4 gnuplot/5.4.0

#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

export GNUTERM=png
for i in {1..435}
do
  ./format_field ../data/july_scans/s01/output_$i.nc 12
  ./format_rate ../data/july_scans/s01/rate_$i.nc
  gnuplot plot_tail.gnu > ../figs/july/s01/slice_$i.svg
done

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
