#!/bin/bash

for i in {1..819}
do
  ./format_field ../data/july_scans/s02/output_$i.nc 12
  ./format_rate ../data/july_scans/s02/rate_$i.nc
  gnuplot plot_tail.gnu > ../figs/july/s02/slice_$i.svg
done
