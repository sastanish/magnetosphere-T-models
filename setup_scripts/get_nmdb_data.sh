#!/bin/bash

PROJ_HOME=$(git rev-parse --show-toplevel)

#for iid in $($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk '/"ID"]/ {print $2}')
for iid in $(cat ../storms.json | jq ".[].ID" | sed 's/"//g')
do

  name=$(cat $PROJ_HOME/storms.json | jq ".$iid.name" | sed 's/"//g')

  syear=$(cat ../storms.json | jq ".$iid.start_date.year" | sed 's/"//g')
  smonth=$(cat ../storms.json | jq ".$iid.start_date.month" | sed 's/"//g')
  sday=$(cat ../storms.json | jq ".$iid.start_date.day" | sed 's/"//g')
  shour=$(cat ../storms.json | jq ".$iid.start_date.hour" | sed 's/"//g')

  eyear=$(cat ../storms.json | jq ".$iid.end_date.year" | sed 's/"//g')
  emonth=$(cat ../storms.json | jq ".$iid.end_date.month" | sed 's/"//g')
  eday=$(cat ../storms.json | jq ".$iid.end_date.day" | sed 's/"//g')
  ehour=$(cat ../storms.json | jq ".$iid.end_date.hour" | sed 's/"//g')

  wget -np -O $PROJ_HOME/data/$name/nmdb/nmdb_data_$name.lst "https://www.nmdb.eu/nest/draw_graph.php?wget=1&allstations&output=ascii&tabchoice=ori&dtype=corr_for_efficiency&date_choice=bydate&start_year=$syear&start_month=$smonth&start_day=$sday&start_hour=$shour&start_min=00&end_year=$eyear&end_month=$emonth&end_day=$eday&end_hour=$ehour&end_min=59&yunits=0"

done
