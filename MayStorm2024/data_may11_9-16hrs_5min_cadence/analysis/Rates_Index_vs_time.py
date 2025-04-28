import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

omni_data = np.genfromtxt("./may11_09h-15h_5min_cadence.dat",dtype=None)

symh = np.zeros(len(omni_data))
rate = np.zeros((len(omni_data),3))

times = []

for i,line in enumerate(omni_data):
    
    time = pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")
    times.append(time)
    time = str(time).replace(" ","_")
    ds_t1 = xr.open_dataset(f"../metrics/2024-05-11-9hrs-16hrs_c2_t1_{i}.nc")
    ds_t2 = xr.open_dataset(f"../metrics/2024-05-11-9hrs-16hrs_c2_t2_{i}.nc")
    ds_t3 = xr.open_dataset(f"../metrics/2024-05-11-9hrs-16hrs_c2_t3_{i}.nc")

    ds_t1 = ds_t1.where(np.sqrt(ds_t1.x**2+ds_t1.y**2+ds_t1.z**2)>1.2)
    ds_t2 = ds_t2.where(np.sqrt(ds_t1.x**2+ds_t1.y**2+ds_t1.z**2)>1.2)
    ds_t3 = ds_t3.where(np.sqrt(ds_t1.x**2+ds_t1.y**2+ds_t1.z**2)>1.2)

    # Determine the maximum reconnection rates
    rate[i,0] = np.abs(ds_t1.c2_t1).max().values
    rate[i,1] = np.abs(ds_t2.c2_t2).max().values
    rate[i,2] = np.abs(ds_t3.c2_t3).max().values

    # Record SymH index
    symh[i] = line[12]

    ds_t1.close()
    ds_t2.close()
    ds_t3.close()

rate_dtimes = pd.DatetimeIndex(times)

# Record kp and DST index
ind_data = np.genfromtxt("./indicies.dat", dtype=None)
kp = np.zeros(len(ind_data))
dst = np.zeros(len(ind_data))
times = []
for i,line in enumerate(ind_data):
    kp[i] = line[3]/10
    dst[i] = line[4]
    time = pd.to_datetime(str(line[0]) + "/" + str(line[1]) + "/" + str(line[2]), format="%Y/%j/%H")
    times.append(time)

ind_times = pd.DatetimeIndex(times)

# Now we graph
fig, (rax1, sax, kax) = plt.subplots(figsize=(10,12),nrows=3)

fig.subplots_adjust(right=(0.8))

rax2 = rax1.twinx()
rax3 = rax1.twinx()
rax3.spines.right.set_position( ("axes", 1.15) )
dax = kax.twinx()

r1plot = rax1.plot(rate_dtimes,rate[:,0],label="t1 rate",color="tab:blue")
r2plot = rax2.plot(rate_dtimes,rate[:,1],label="t2 rate",color="tab:red")
r3plot = rax3.plot(rate_dtimes,rate[:,2],label="t3 rate",color="tab:green")

splot = sax.plot(rate_dtimes,symh,label="Symh",color="tab:red")

kplot = kax.plot(ind_times,kp,label="Kp",color="tab:green")
dplot = dax.plot(ind_times,dst,label="DST",color="tab:blue")

rax1.set(ylabel="t1")
rax2.set(ylabel="t2")
rax3.set(ylabel="t3")

sax.set(ylabel="Symh")
dax.set(ylabel="DST")

kax.set(ylabel="Kp")

fig.legend()
plt.savefig("Rates_Index_vs_time.png")
