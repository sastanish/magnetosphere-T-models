#!/bin/bash

declare -a years=(2018 2022 2023 2024 2024 2021 2022 2023 2024 2024 2022 2023 2023 2024)
declare -a months=(08 03 03 03 08 11 10 04 03 10 01 02 11 05)
declare -a days=(25 13 22 03 11 03 22 23 24 10 14 26 06 10)
# loop through all dates
for i in {0..13};
do
  wget -np -O "nmdb_data_${years[$i]}-${months[$i]}-${days[$i]}.lst" "https://www.nmdb.eu/nest/draw_graph.php?wget=1&allstations&output=ascii&tabchoice=ori&dtype=corr_for_efficiency&date_choice=bydate&start_year=${years[$i]}&start_month=${months[$i]}&start_day=$((${days[$i]}-1))&start_hour=00&start_min=00&end_year=${years[$i]}&end_month=${months[$i]}&end_day=$((${days[$i]}+2))&end_hour=23&end_min=59&yunits=0"
done
