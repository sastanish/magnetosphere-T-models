#!/bin/bash

for i in {1..435}
do
  ./format_field ../data/july_scans/s01/output_$i.nc 12
  ./format_rate ../data/july_scans/s01/rate_$i.nc
  gnuplot plot_tail.gnu > ../figs/july/s01/slice_$i.svg
done
