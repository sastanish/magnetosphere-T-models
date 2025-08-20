#!/bin/bash

# loop through all dates
mkdir figs
cd figs
for iid in $(../../input_data/json.sh < ../../input_data/storms.json | awk '/"ID"]/ {print $2}')
do
  idnq=$(echo $iid | sed 's/"//g')
  mkdir $idnq
done
