#!/bin/bash

# loop through all dates
for iid in $(../input_data/json.sh < ../input_data/storms.json | awk '/"ID"]/ {print $2}')
do
  idnq=$(echo $iid | sed 's/"//g')
  rm -r $idnq
  mkdir $idnq
  cp ../input_data/data/"$idnq"_TA16_parameters.lst $idnq/input_data.lst
  cp ./TA16_RBF.par $idnq/TA16_RBF.par
  cp ./input_parameters.txt $idnq/input_parameters.txt
  cp ./ta16.out $idnq/ta16.out
  cp ./recon.out $idnq/recon.out
done
