#!/bin/bash

# loop through all dates
for date in 2018-08-25 2022-03-13 2023-03-22 2024-03-03 2024-08-11 2021-11-03 2022-10-22 2023-04-23 2024-03-24 2024-10-10 2022-01-14 2023-02-26 2023-11-06 2024-05-10
do
  echo "
<details>
<summary>Storm ${date}</summary>
  <figure>
    <img src=\"../figs/TA16/${date}/x-point_location.png\">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src=\"../figs/TA16/${date}/protons_neutrons_and_x-points.png\">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
</details>" >> hold.md
done
