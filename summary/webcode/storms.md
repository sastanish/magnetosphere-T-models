---
title: Near Earth x-point
date: 2025/06/11
author: Sage Stanish
css: theme.css
---

This document details the solar storms that we have simulated using the TA16 empirical magnetosphere of Tsyganenko. 

# A note to David

All of these storms have been re-run using the correct index. Our N-index lines up well with Tsyganenko's with some deviation, perhaps due to the cadence of our underlying data.

<figure>
  <img src="../figs/my_inds_vs_Tsyganenko_2018.png">
  <figcaption>A comparison between Tsyganenko's (blue) and my (red) input data: $\langle Nindex \rangle$, $\langle SymHc \rangle$, and $\langle BZ \rangle$.</figcaption>
</figure>

# Overview of storms

The following storms are all simulated over the domain $x,y,z \in (-17,0) \times (-10,10) \times (-4,4)$ with a resolution of 5cells per earth radii in $x,y$ and 10 cells per radii in $z$. 
<details>
<summary>Storm 2018-08-25</summary>
  ## Notes
  - Could be a potential x-point that forms within 6Re, but need to look in higher resolution.
  <figure>
    <img src="../figs/TA16/2018-08-25/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2018-08-25/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2018-08-25/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
  <figure>
    <img src="../figs/TA16/2018-08-25/potential_x-point.png">
    <figcaption>A potential x-point forming during the 2018 storm.</figcaption>
  </figure>
</details>

<details>
<summary>Storm 2021-11-03</summary>
  ## Notes
  - The current sheet gets quite close to the earth, but there is no x-point formation.
  <figure>
    <img src="../figs/TA16/2021-11-03/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2021-11-03/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2021-11-03/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2022-01-14</summary>
  <figure>
    <img src="../figs/TA16/2022-01-14/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2022-01-14/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2022-01-14/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2022-03-13</summary>
  <figure>
    <img src="../figs/TA16/2022-03-13/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2022-03-13/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2022-03-13/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2022-10-22</summary>
  <figure>
    <img src="../figs/TA16/2022-10-22/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2022-10-22/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2022-10-22/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2023-02-26</summary>
  <figure>
    <img src="../figs/TA16/2023-02-26/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2023-02-26/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2023-02-26/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2023-03-22</summary>
  <figure>
    <img src="../figs/TA16/2023-03-22/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2023-03-22/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2023-03-22/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2023-04-23</summary>
  <figure>
    <img src="../figs/TA16/2023-04-23/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2023-04-23/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2023-04-23/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2023-11-06</summary>
  <figure>
    <img src="../figs/TA16/2023-11-06/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2023-11-06/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2023-11-06/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2024-03-24</summary>
  <figure>
    <img src="../figs/TA16/2024-03-24/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2024-03-24/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2024-03-24/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2024-03-03</summary>
  <figure>
    <img src="../figs/TA16/2024-03-03/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2024-03-03/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2024-03-03/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2024-05-10</summary>
  <figure>
    <img src="../figs/TA16/2024-05-10/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2024-05-10/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2024-05-10/x-point_tail_from_above_MayStorm.png">
    <img src="../figs/TA16/2024-05-10/x-point_tail_MayStorm.png">
    <figcaption>These images show the May Storm during at the nearest point of the x-point. The colored surface shows a contour of the reconnection rate at ~75 while the streamlines trace the field. The x-point extends all around the night-side of the earth.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2024-05-10/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2024-08-11</summary>
  <figure>
    <img src="../figs/TA16/2024-08-11/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2024-08-11/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2024-08-11/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

<details>
<summary>Storm 2024-10-10</summary>
  <figure>
    <img src="../figs/TA16/2024-10-10/x-point_location.png">
    <figcaption>X-point location as determined by the pressure (red) and critical (blue) methods. This is compared with some input parameters in SymHc, N-Index, and BZ (averaged in red).</figcaption>
  </figure>
  <figure>
    <img src="../figs/TA16/2024-10-10/protons_neutrons_and_x-points.png">
    <figcaption>X-point location in comparison with neutron monitor data and GOES proton flux. The neutron monitors have been normalized to the unit interval. The vertical line corresponds to the peak in the proton flux.</figcaption>
  </figure>
  <p>
  This movie shows a slice through the magneto-tail of the earth at $y=0$. The background shows the magnitude of the reconnection rate with the magnetic field projected onto it in white.
  </p>
  <video controls>
    <source src="../figs/TA16/2024-10-10/tail_slices/tail_movie.mp4" type=video/mp4>
  </video>
</details>

# Methods

To determine the location of the x-point in these simulations, we need to identify dips in the pressure that correspond to the x-point and surrounding tail-current. There are two methods of determining the location of the x-point that we use in data below. The first, in red, is the pressure ball method and the second, in blue, is the critical point method. 

## Pressure-Ball Method.
The algorithm is as follows:

 1) Find the global minimum in $P_B=|B|^2$.
 2) Let $\mathcal{B}(P_{B,\min})$ be the ball of points of radius $0.5Re$ around the minima of the magnetic pressure.
 3) Find the maximum reconnection rate within $\mathcal{B}(P_{B,\min})$. This rate corresponds to either the center of the x-point or the center of the current sheet/flux rope. Record this rate as $R_i$ and location as $x_i$.
 4) Repeat this process 3 times, removing $\mathcal{B}(P_{B,\min})$ from the pressure each time. Then the reconnection rate that we keep is the one that is closest to earth, $\min(|x_i|)$. 

This repeated procedure pics out the closest drop in pressure to the earth. We iterate this process since the current sheet will often have a lower magnetic pressure than the x-point itself. So, we have to ensure that we get the actual x-point by considering the closest such pressure drop.

Note: This method will always return some point, even if no x-point or current sheet exists.

## Critical Point Method
The algorithm:

 1) Consider a North-South slice along the magneto-tail ($y=0$)
 2) Via rolling array operations, find all points (not within the Earth) where both derivatives of the pressure ($\partial_x P=\partial_z P=0$) change sign. 
 3) Consider the closest point the center of the x-point and record it's location and reconnection rate.

This method does not consider the full 3d-field, but does reproduce the x-point location well. It agrees with the Pressure-Ball method when an x-point is present and also returns *NaN* for times when no x-point exists.


