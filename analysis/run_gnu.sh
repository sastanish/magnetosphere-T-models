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
#SBATCH --time=01:00:00
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

module load netcdf-fortran/intel-2020.4 intel/intel-2020.4

#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

export GNUTERM=dumb
for i in {1..699}
do
  ./format_field.out ../data/july_scans/s06/output_$i.nc 12
  ./format_rate.out ../data/july_scans/s06/rate_$i.nc
  gnuplot plot_tail.gnu > ../figs/july/s06/slice_$i.eps
done

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
