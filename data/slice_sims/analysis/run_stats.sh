#!/bin/bash



export OMP_NUM_THREADS=8

# loop through all dates
for iid in $(../../input_data/json.sh < ../../input_data/storms.json | awk '/"ID"]/ {print $2}')
do
  idnq=$(echo $iid | sed 's/"//g')
  nlines=$(cat ../${idnq}/input_data.lst | wc -l)
  echo ${nlines}
  ./x-point_location.out ../$idnq/ 1 ${nlines}
  ./field_drop.out ../$idnq/ 1 ${nlines}
  ./flux.out ../$idnq/ 1 ${nlines}
  ./pressure.out ../$idnq/ 1 ${nlines}
done
