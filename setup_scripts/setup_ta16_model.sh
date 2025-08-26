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
    cp "$PROJ_HOME/data/$name/TS05"
  fi 

done
