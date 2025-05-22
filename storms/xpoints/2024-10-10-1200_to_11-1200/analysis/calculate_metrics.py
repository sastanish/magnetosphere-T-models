import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

omni_data = np.genfromtxt("../input_data.lst",dtype=None)
Nt = len(omni_data)

symh = np.zeros(Nt)
rate = np.zeros(Nt)
rad = np.zeros(Nt)
t1_rate = np.zeros(Nt)
t2_rate = np.zeros(Nt)
t3_rate = np.zeros(Nt)
t1_rad = np.zeros(Nt)
t2_rad = np.zeros(Nt)
t3_rad = np.zeros(Nt)
times = []

for i,line in enumerate(omni_data):
    
    times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

    ds = xr.open_dataset(f"../data/2024_data_{i+1}.nc")
    ds = ds.where(ds.z<=-2)

    totalrate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3
    totalrate = totalrate.where(totalrate==totalrate.max(), drop=True).squeeze()
    rate[i] = totalrate.values
    rad[i] = np.sqrt(totalrate.x**2+totalrate.y**2+totalrate.z**2).values

    t_rate = ds.mag_c2_t1.where(ds.mag_c2_t1==ds.mag_c2_t1.max(), drop=True).squeeze()
    t1_rate[i] = t_rate.values
    t1_rad[i] = np.sqrt(t_rate.x**2+t_rate.y**2+t_rate.z**2).values

    t_rate = ds.mag_c2_t2.where(ds.mag_c2_t2==ds.mag_c2_t2.max(), drop=True).squeeze()
    t2_rate[i] = t_rate.values
    t2_rad[i] = np.sqrt(t_rate.x**2+t_rate.y**2+t_rate.z**2).values

    t_rate = ds.mag_c2_t3.where(ds.mag_c2_t3==ds.mag_c2_t3.max(), drop=True).squeeze()
    t3_rate[i] = t_rate.values
    t3_rad[i] = np.sqrt(t_rate.x**2+t_rate.y**2+t_rate.z**2).values

    ds.close()

np.savetxt('totalRate.txt', rate, header='Maximum of the total reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('totalRad.txt', rad, header='Radius of the maximum total reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('t1Rate.txt', t1_rate, header='Maximum of first term in the reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('t1Rad.txt', t1_rad, header='Radius of the maximum first term in the reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('t2Rate.txt', t2_rate, header='Maximum of second term in the reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('t2Rad.txt', t2_rad, header='Radius of the maximum second term in the reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('t3Rate.txt', t3_rate, header='Maximum of third term in the reconnection rate during the 2023/10/10 - 2023/10/11 storm')
np.savetxt('t3Rad.txt', t3_rad, header='Radius of the maximum third term in the reconnection rate during the 2023/10/10 - 2023/10/11 storm')
