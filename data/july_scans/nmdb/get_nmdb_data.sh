#!/bin/bash

for iid in $(../input_data/json.sh < ../input_data/storms.json | awk '/"ID"]/ {print $2}')
do
  # Parse the start time
  syear=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"start_date","year"/ { print $2 }'| sed 's/"//g')
  smonth=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"start_date","month"/ { print $2 }'| sed 's/"//g')
  sday=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"start_date","day"/ { print $2 }'| sed 's/"//g')
  shour=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"start_date","hour"/ { print $2 }'| sed 's/"//g')

  # Now the end time
  eyear=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"end_date","year"/ { print $2 }'| sed 's/"//g')
  emonth=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"end_date","month"/ { print $2 }'| sed 's/"//g')
  eday=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"end_date","day"/ { print $2 }'| sed 's/"//g')
  ehour=$(../input_data/json.sh < ../input_data/storms.json | awk '/'$iid',"end_date","hour"/ { print $2 }'| sed 's/"//g')

  echo $shour
  echo $ehour

  # Now we curl in the data
  idnq=$(echo $iid | sed 's/"//g')
  wget -np -O nmdb_data_$idnq.lst "https://www.nmdb.eu/nest/draw_graph.php?wget=1&allstations&output=ascii&tabchoice=ori&dtype=corr_for_efficiency&date_choice=bydate&start_year=$syear&start_month=$smonth&start_day=$sday&start_hour=$shour&start_min=00&end_year=$eyear&end_month=$emonth&end_day=$eday&end_hour=$ehour&end_min=59&yunits=0"
done
