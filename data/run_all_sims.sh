#!/bin/bash


PROJ_HOME=$(git rev-parse --show-toplevel)

export FOR_COARRAY_NUM_IMAGES=15
for iid in $($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk '/"ID"]/ {print $2}')
do

  # Parse the short name
  #name=$(echo $iid | sed 's/"//g')
  #name=$($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk "/${iid},name/ {print \"$0\"}")

  include=$(cat $PROJ_HOME/storms.json | jq ".$iid.include" | sed 's/"//g')
  name=$(cat $PROJ_HOME/storms.json | jq ".$iid.name" | sed 's/"//g')

  if [ $include = 'yes' ]; then
    cd $PROJ_HOME/data/$name/TA16
    ./pressure_and_flux.out
  fi 
  cd $PROJ_HOME/data
done


export OMP_NUM_THREADS=15
cd $PROJ_HOME/data/Feb2022/TA16
nlines=$(wc -l < input_data.lst)
./main_TA16.out 1 $nlines
./recon.out 1 $nlines
cd $PROJ_HOME/data/Jun2015/TA16
nlines=$(wc -l < input_data.lst)
./main_TA16.out 1 $nlines
./recon.out 1 $nlines
