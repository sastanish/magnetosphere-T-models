#!/bin/bash

cd ../figs
# loop through all dates
for iid in $(../../input_data/json.sh < ../../input_data/storms.json | awk '/"ID"]/ {print $2}')
do
  idnq=$(echo $iid | sed 's/"//g')
  cd $idnq/slices
  ffmpeg -i tail_slice_%d.png -framerate 12 ../${idnq}_movie.mp4
  cd ../../
done
