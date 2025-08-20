#!/bin/bash

export OMP_NUM_THREADS=10

# loop through all dates
for iid in $(../input_data/json.sh < ../input_data/storms.json | awk '/"ID"]/ {print $2}')
do
  idnq=$(echo $iid | sed 's/"//g')
  cd $idnq
  nlines=$(cat ./input_data.lst | wc -l)
  echo ${nlines}
  ./ta16.out 1 ${nlines}
  ./recon.out 1 ${nlines}
  cd ../
done
