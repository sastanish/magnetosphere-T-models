#!/bin/zsh

for date in 2018-08-25_to_27 2023-03-23_to_24 2021-11-03-1800_to_04-1800 2024-10-10-1200_to_11-1200 2022-03-13-10_to_14-08
do
  # setup data directory and move input data in
  mkdir $date
  cp storm.py $date/storm.py
  cp TA16_RBF.par $date/TA16_RBF.par
  cd $date
  mkdir data
  cp /home/us/phd/code/magnetosphere-T-models/data_handling/data/omni_?min_($date)_TA16_parameters.lst input_data.lst
  cd ../
done
