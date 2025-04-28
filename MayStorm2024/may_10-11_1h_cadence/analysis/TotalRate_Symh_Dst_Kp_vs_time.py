import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

omni_data = np.genfromtxt("../1hour_cadence.dat",dtype=None)

symh = np.zeros(len(omni_data))
rate = np.zeros((len(omni_data),3))

times = []

for i,line in enumerate(omni_data):
    
    time = pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")
    times.append(time)
    time = str(time).replace(" ","_")
    ds_t1 = xr.open_dataset("../data/" + str(time) + "_c2_t1.nc")
    ds_t2 = xr.open_dataset("../data/" + str(time) + "_c2_t2.nc")
    ds_t3 = xr.open_dataset("../data/" + str(time) + "_c2_t3.nc")

    ds_t1 = ds_t1.where(np.sqrt(ds_t1.x**2+ds_t1.y**2+ds_t1.z**2)>1.2)
    ds_t2 = ds_t2.where(np.sqrt(ds_t1.x**2+ds_t1.y**2+ds_t1.z**2)>1.2)
    ds_t3 = ds_t3.where(np.sqrt(ds_t1.x**2+ds_t1.y**2+ds_t1.z**2)>1.2)

    # Determine the maximum reconnection rates
    rate[i,0] = np.abs(ds_t1.c2_t1.sel(x=slice(-14.8,1.8),y=slice(-3.9,3.9),z=slice(-3.9,3.9))).max().values
    rate[i,1] = np.abs(ds_t2.c2_t2.sel(x=slice(-14.8,1.8),y=slice(-3.9,3.9),z=slice(-3.9,3.9))).max().values
    rate[i,2] = np.abs(ds_t3.c2_t3.sel(x=slice(-14.8,1.8),y=slice(-3.9,3.9),z=slice(-3.9,3.9))).max().values

    # Record SymH index
    symh[i] = line[10]

    ds_t1.close()
    ds_t2.close()
    ds_t3.close()

rate_dtimes = pd.DatetimeIndex(times)

# Record kp index
kp_data = np.genfromtxt("./kp_may_2024.dat", dtype=None)
kp = np.zeros(len(kp_data))
times = []
for i,line in enumerate(kp_data):
    kp[i] = line[7]
    time = pd.to_datetime(str(line[0]) + "/" + str(line[1]) + "/" + str(line[2]), yearfirst=False)
    times.append(time)

kp_times = pd.DatetimeIndex(times)

# Record DST index
dst_data = np.genfromtxt("./dst_may_2024.dat", dtype=None, skip_header=18)
dst = np.zeros(len(dst_data))
times = []
for i,line in enumerate(dst_data):
    dst[i] = line[3]
    time = pd.to_datetime(str(line[0]) + "_" + str(line[1]), format="%Y-%m-%d_%H:%M:%S.000")
    times.append(time)

dst_times = pd.DatetimeIndex(times)

# Now we graph
fig, rax = plt.subplots(figsize=(10,4))

fig.subplots_adjust(right=(0.75))

sax = rax.twinx()
dax = rax.twinx()
dax.spines.right.set_position( ("axes", 1.1) )
kax = rax.twinx()
kax.spines.right.set_position( ("axes", 1.23) )

rplot = rax.plot(rate_dtimes,np.sum(rate,axis=-1),label="rate",color="tab:blue")
splot = sax.plot(rate_dtimes,symh,label="Symh",color="tab:red")
dplot = dax.plot(dst_times,dst,label="DST",color="tab:orange")
kplot = kax.plot(kp_times,kp,label="Kp",color="tab:green")

rax.set(ylabel="rate")
sax.set(ylabel="Symh")
dax.set(ylabel="DST")
kax.set(ylabel="Kp")

fig.legend()
plt.savefig("TotalRate_Symh_Dst_Kp_vs_time.png")

