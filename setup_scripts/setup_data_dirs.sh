#!/bin/bash


PROJ_HOME=$(git rev-parse --show-toplevel)

for iid in $($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk '/"ID"]/ {print $2}')
do

  # Parse the short name
  #name=$(echo $iid | sed 's/"//g')
  #name=$($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk "/${iid},name/ {print \"$0\"}")

  include=$(cat $PROJ_HOME/storms.json | jq ".$iid.include" | sed 's/"//g')
  name=$(cat $PROJ_HOME/storms.json | jq ".$iid.name" | sed 's/"//g')

  if [ $include = 'yes' ]; then
    mkdir -p "$PROJ_HOME/data/$name/TS05"
    mkdir -p "$PROJ_HOME/data/$name/TA16"
    mkdir -p "$PROJ_HOME/data/$name/omni"
    mkdir -p "$PROJ_HOME/data/$name/nmdb"
    cp $PROJ_HOME/sim_codes/main_TA16.out $PROJ_HOME/data/$name/TA16/
    cp $PROJ_HOME/sim_codes/recon.out $PROJ_HOME/data/$name/TA16/
    cp $PROJ_HOME/sim_codes/pressure_and_flux.out $PROJ_HOME/data/$name/TA16/
    cp $PROJ_HOME/sim_codes/TA16_RBF.par $PROJ_HOME/data/$name/TA16/
    cp $PROJ_HOME/sim_codes/input_parameters.txt $PROJ_HOME/data/$name/TA16/
  fi 

done
