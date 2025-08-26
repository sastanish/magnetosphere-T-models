#!/bin/bash

PROJ_HOME=$(git rev-parse --show-toplevel)

#for iid in $($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk '/"ID"]/ {print $2}')
for iid in $(cat ../../storms.json | jq ".[].ID" | sed 's/"//g')
do

  name=$(cat $PROJ_HOME/storms.json | jq ".$iid.name" | sed 's/"//g')
  
  cp data/omni_data_$name* $PROJ_HOME/data/$name/omni
  cp data/${name}_TA16_parameters.lst $PROJ_HOME/data/$name/omni

done
