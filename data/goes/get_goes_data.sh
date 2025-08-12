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

  stime="${syear}-${smonth}-${sday}T${shour}:00:00"
  echo $stime
  etime="${eyear}-${emonth}-${eday}T${ehour}:59:00"
  echo $etime

  # Now we curl in the data
  idnq=$(echo $iid | sed 's/"//g')
  wget -np -O goes_flux_data_$idnq.lst "https://lasp.colorado.edu/space-weather-portal/latis/dap/goesp_part_flux_P5M.asc?time>=${stime}&time<=${etime}&formatTime(yyyy-MM-dd HH:mm:ss)"
done
