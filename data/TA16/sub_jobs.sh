#!/bin/bash

# loop through all dates
for date in 2018-08-25 2022-03-13 2023-03-22 2024-03-03 2021-11-03 2022-10-22 2023-04-23 2022-01-14 2023-02-26 2023-11-06
do
  cd $date
  sbatch ./run_storm.sh
  cd ../
done
