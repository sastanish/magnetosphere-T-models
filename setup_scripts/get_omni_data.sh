#!/bin/bash

PROJ_HOME=$(git rev-parse --show-toplevel)

#for iid in $($PROJ_HOME/setup_scripts/json.sh < $PROJ_HOME/storms.json | awk '/"ID"]/ {print $2}')
for iid in $(cat ../storms.json | jq ".[].ID" | sed 's/"//g')
do

  name=$(cat $PROJ_HOME/storms.json | jq ".$iid.name" | sed 's/"//g')

  year=$(cat ../storms.json | jq ".$iid.start_date.year" | sed 's/"//g')
  month=$(cat ../storms.json | jq ".$iid.start_date.month" | sed 's/"//g')
  day=$(cat ../storms.json | jq ".$iid.start_date.day" | sed 's/"//g')
  hour=$(cat ../storms.json | jq ".$iid.start_date.hour" | sed 's/"//g')
  start_date=$(echo "$year$month$day$hour")

  year=$(cat ../storms.json | jq ".$iid.end_date.year" | sed 's/"//g')
  month=$(cat ../storms.json | jq ".$iid.end_date.month" | sed 's/"//g')
  day=$(cat ../storms.json | jq ".$iid.end_date.day" | sed 's/"//g')
  hour=$(cat ../storms.json | jq ".$iid.end_date.hour" | sed 's/"//g')
  end_date=$(echo "$year$month$day$hour")

  curl -d "activity=retrieve&res=5min&spacecraft=omni_5min&start_date=$start_date&end_date=$end_date&vars=14&vars=15&vars=16&vars=22&vars=23&vars=24&vars=25&vars=26&vars=27&vars=37&vars=41&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi | sed 1,18d | tac | sed 1,15d | tac > $PROJ_HOME/setup_scripts/format_omni_data/data/omni_data_$name.lst

done
