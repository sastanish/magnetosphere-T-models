#!/bin/bash

# loop through all dates
for date in 2024-05-10
# 2018-08-25 2022-03-13 2023-03-22 2024-03-03 2024-08-11 2021-11-03 2022-10-22 2023-04-23 2024-03-24 2024-10-10 2022-01-14 2023-02-26 2023-11-06 2024-05-10
do
  mkdir $date
  cp ../prepared_omni_data/TA16/data/${date}_TA16_parameters.lst $date/input_data.lst
  cp ./TA16_RBF.par $date/TA16_RBF.par
  cp ./storm.py $date/storm.py
  cd $date
  ln -s ../../../sim_codes/src/TsyganenkoWrapper.cpython-312-x86_64-linux-gnu.so .
  cd ../
done
