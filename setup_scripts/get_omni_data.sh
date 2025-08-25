#!/bin/bash

for iid in $(./json.sh < storms.json | awk '/"ID"]/ {print $2}')
do
  # Parse the start time
  year=$(./json.sh < storms.json | awk '/'$iid',"start_date","year"/ { print $2 }')
  month=$(./json.sh < storms.json | awk '/'$iid',"start_date","month"/ { print $2 }')
  day=$(./json.sh < storms.json | awk '/'$iid',"start_date","day"/ { print $2 }')
  hour=$(./json.sh < storms.json | awk '/'$iid',"start_date","hour"/ { print $2 }')
  start_date=$(echo "$year$month$day$hour" | sed 's/"//g')

  # Now the start time
  year=$(./json.sh < storms.json | awk '/'$iid',"end_date","year"/ { print $2 }')
  month=$(./json.sh < storms.json | awk '/'$iid',"end_date","month"/ { print $2 }')
  day=$(./json.sh < storms.json | awk '/'$iid',"end_date","day"/ { print $2 }')
  hour=$(./json.sh < storms.json | awk '/'$iid',"end_date","hour"/ { print $2 }')
  end_date=$(echo "$year$month$day$hour" | sed 's/"//g')
  echo $end_date

  # Now we curl in the data
  idnq=$(echo $iid | sed 's/"//g')
  curl -d "activity=retrieve&res=5min&spacecraft=omni_5min&start_date=$start_date&end_date=$end_date&vars=14&vars=15&vars=16&vars=22&vars=23&vars=24&vars=25&vars=26&vars=27&vars=37&vars=41&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi | sed 1,18d | tac | sed 1,15d | tac > ./data/omni_data_$idnq.lst

done
