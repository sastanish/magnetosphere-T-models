# Main Takeaways

 1) Only 4 out of the 13 storms studied show near earth reconnection
 2) Only 2 out of these 4 show an enhancement in neutron showers from detectors corresponding to the x-point formation
 3) Neutron enhancement occurs only when the reconnection region is within 5 earth radii. Regions farther from the earth have their particles slowed down before they impact the atmosphere.
 4) All of the reconnection rates are the same order at $\approx \eta 10^{2}$.


## Storm List

The following is a list of solar storms that we have simulated using one of Tsyganenko's magnetosphere models. For each storm, I have listed the simulated interval, cadence, resolution, model, and qualitative state of the magneto-tail. 

| Dates | Cadence | Domain | Resolution | Model | Tail-State |
| ---   | ---     |    --- |        --- | ---   |        --- |
| 2018/08/25 00:30 - 27 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2021/11/03 18:30 - 04 18:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2022/01/14 01:00 - 16 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2022/10/22 00:30 - 12 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2023/02/26 12:30 - 27 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2023/03/23 00:30 - 24 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2023/04/24 01:20 - 24 07:59 | 1min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/30$ | TA16 | Xpoint |
| 2023/11/05 00:30 - 06 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | Disturbed Tail |
| 2024/03/03 00:30 - 04 12:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | No X-point |
| 2024/03/24 14:31 - 24 18:00 | 1min | $x,y,z \in (-8,0) \times (-1,1) \times (0,2)$   | $1/30$ | TA16 | Xpoint |
| 2024/05/10 18:00 - 06 22:40 | 1min | $x,y,z \in (-7,0) \times (-1,1) \times (0,2)$   | $1/30$ | TA16, TS05 | Xpoint |
| 2024/08/11 00:30 - 13 23:40 | 10min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$ | $1/20$ | TA16 | Disturbed Tail |
| 2024/10/10 15:00 - 11 10:59 | 1min | $x,y,z \in (-8,0) \times (-2,2) \times (-2,2)$   | $1/25$ | TA16, TS05 | Xpoint |

## Storms with Near Earth X-Points

### Methods
To determine the location of the x-point in these simulations, we need to identify dips in the pressure that correspond to the x-point and surrounding tail-current. There are two methods of determining the location of the x-point that we use in data below. The first, in red, is the pressure ball method and the second, in blue, is the critical point method. 

#### Pressure-Ball Method.
The algorithm is as follows:

 1) Find the global minimum in $P_B=|B|^2$.
 2) Let $\mathcal{B}(P_{B,\min})$ be the ball of points of radius $0.5Re$ around the minima of the magnetic pressure.
 3) Find the maximum reconnection rate within $\mathcal{B}(P_{B,\min})$. This rate corresponds to either the center of the x-point or the center of the current sheet/flux rope. Record this rate as $R_i$ and location as $x_i$.
 4) Repeat this process 3 times, removing $\mathcal{B}(P_{B,\min})$ from the pressure each time. Then the reconnection rate that we keep is the one that is closest to earth, $\min(|x_i|)$. 

This repeated procedure pics out the closest drop in pressure to the earth. We iterate this process since the current sheet will often have a lower magnetic pressure than the x-point itself. So, we have to ensure that we get the actual x-point by considering the closest such pressure drop.

Note: This method will always return some point, even if no x-point or current sheet exists.

#### Critical Point Method
The algorithm:

 1) Consider a North-South slice along the magneto-tail ($y=0$)
 2) Via rolling array operations, find all points (not within the Earth) where both derivatives of the pressure ($\partial_x P=\partial_z P=0$) change sign. 
 3) Consider the closest point the center of the x-point and record it's location and reconnection rate.

This method does not consider the full 3d-field, but does reproduce the x-point location well. It agrees with the Pressure-Ball method when an x-point is present and also returns *NaN* for times when no x-point exists.

### April 2023
The following is data from the April 24, 2023 storm. This storm forms a clear x-point between 01:00 - 08:00. 

![](./2023-04-23-1200_to_24-2400/analysis/dip_vs_neutrons.png)

#### Features:

 - An x-point forms from current sheet incursion. It is steady, but never passes within 5 earth radii.
 - The actual reconnection event is short around the closest point.
 - Little to no Neutron enhancement. There may be some enhancement around the spike in reconnection around 2am, however this could be noise. 

#### Images:

![01:20](./2023-04-23-1200_to_24-2400/time_0120.png)
![01:45](./2023-04-23-1200_to_24-2400/time_0145.png)
![02:05](./2023-04-23-1200_to_24-2400/time_0205.png)
![02:30](./2023-04-23-1200_to_24-2400/time_0230.png)
![04:12](./2023-04-23-1200_to_24-2400/time_0412.png)
![07:05](./2023-04-23-1200_to_24-2400/time_0705.png)

### March 2024

![](./2024-03-24_to_25/analysis/dip_vs_neutrons.png)

### May 2024
The following is data from the May 10-11, 2024 storms. These storms forms a clear x-point between 01:00 - 08:00. 

![](./2024-05-10_to_11/analysis/dip_vs_neutrons.png)

#### Features:

 - Strong x-point forms that oscillates between a strong, close point before being rotated out of view northward.
 - Neutron enhancement during closest x-point times.

#### Images:

![10-18:00](./2024-05-10_to_11/time_10_1800.png)
![10-18:20](./2024-05-10_to_11/time_10_1820.png)
![10-23:00](./2024-05-10_to_11/time_10_2300.png)
![11-00:00](./2024-05-10_to_11/time_11_0000.png)
![11-01:00](./2024-05-10_to_11/time_11_0100.png)
![11-01:45](./2024-05-10_to_11/time_11_0145.png)
![11-03:00](./2024-05-10_to_11/time_11_0300.png)
![11-07:20](./2024-05-10_to_11/time_11_0720.png)
![11-07:30](./2024-05-10_to_11/time_11_0730.png)
![11-12:00](./2024-05-10_to_11/time_11_1200.png)


### October 2024

![](./2024-10-10-1200_to_11-1200/analysis/dip_vs_neutrons.png)

#### Images:

![16:04](./2024-10-10-1200_to_11-1200/time_1604.png)
![17:00](./2024-10-10-1200_to_11-1200/time_1700.png)
![18:00](./2024-10-10-1200_to_11-1200/time_1800.png)
![22:00](./2024-10-10-1200_to_11-1200/time_2200.png)
![23:00](./2024-10-10-1200_to_11-1200/time_2300.png)
![01:00](./2024-10-10-1200_to_11-1200/time_0100.png)
![06:00](./2024-10-10-1200_to_11-1200/time_0600.png)
