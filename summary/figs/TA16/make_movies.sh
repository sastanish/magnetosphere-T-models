#!/bin/bash

# loop through all dates
for date in 2018-08-25 2022-03-13 2023-03-22 2024-03-03 2024-08-11 2021-11-03 2022-10-22 2023-04-23 2024-03-24 2024-10-10 2022-01-14 2023-02-26 2023-11-06 2024-05-10
do
  cd $date/tail_slices
  rm tail_movie.mp4
  cat $(find . -name "*.png" -print | sort -V) | ffmpeg -framerate 7 -f image2pipe -i - tail_movie.mp4
  cd ../../
done
