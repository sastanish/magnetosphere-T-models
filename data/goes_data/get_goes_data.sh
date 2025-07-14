#!/bin/bash

# loop through all dates
for date in 2022-03-13 2023-03-22 2024-03-03 2024-08-11 2021-11-03 2022-10-22 2023-04-23 2024-03-24 2024-10-10 2022-01-14 2023-02-26 2023-11-06 2024-05-10
do
  wget -np -O goes_flux_data_${date}.lst "https://lasp.colorado.edu/space-weather-portal/latis/dap/goesp_part_flux_P5M.asc?time>=${date}&limit(864)&formatTime(yyyy-MM-dd HH:mm:ss)"
done
