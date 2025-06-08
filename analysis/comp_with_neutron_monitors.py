import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

# Clean up neutron monitor data
neutron_data = np.genfromtxt("./neutron_data.txt",
                             skip_header=42,dtype=None,names=True,missing_values="null"
                             )
times = []
for i in range(len(neutron_data)):
    times.append(pd.to_datetime(str(neutron_data["DATE"][i]) + "_" + str(neutron_data["HOUR"][i]), format="%Y-%m-%d_%H:%M:%S"))

ntime = pd.to_datetime(times)

omni_data = np.genfromtxt("../input_data.lst",dtype=None)

Nt = len(omni_data)
symh  = np.zeros(Nt)
AeInd = np.zeros(Nt)
BZ = np.zeros(Nt)
times = []

for i in range(Nt):
    line = omni_data[i]
    
    times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))
    symh[i] = line[14]
    AeInd[i] = line[13]
    BZ[i] = line[6]

xtime =pd.to_datetime(times)

times = []
brate = []
brad = []
crate = []
crad = []
for line in np.genfromtxt('x-point_rate.txt',dtype=None):
    times.append(str(line[0]).replace("_"," "))
    brate.append(float(line[1]))
    brad.append( float(line[2]))
    crate.append(float(line[3]))
    crad.append( float(line[4]))

rtime =pd.to_datetime(times)

sel = slice(400,1100)

fig, ax = plt.subplots(nrows=8,figsize=(8,24),sharex=True)

colors = matplotlib.color_sequences["tab20"]
names = neutron_data.dtype.names

for i in range(12):
    ax[0].plot(ntime[4:],moving_average(norm(neutron_data[names[i+2]]),5),label=names[i+2],color=colors[i],alpha=0.9)

ax[0].legend()
ax[0].set_ylabel('monitor data shape')

ax[1].plot(rtime,brate,color="tab:red",linestyle="",marker=".",label="max rate",alpha=0.9)
ax[1].plot(rtime,crate,color="tab:blue",linestyle="",marker=".",label="max rate",alpha=0.9)
ax[1].set_ylabel('rate')
ax[1].set_yscale('log')

for i in range(12):
    ax[2].plot(ntime[4:],moving_average(norm(neutron_data[names[i+14]]),5),label=names[i+14],color=colors[i],alpha=0.9)

ax[2].legend()
ax[2].set_ylabel('monitor data shape')

ax[3].plot(rtime,brad,color="tab:red",linestyle="",marker=".",alpha=0.9)
ax[3].plot(rtime,crad,color="tab:blue",linestyle="",marker=".",alpha=0.9)
ax[3].set_ylabel('radius')

for i in range(10):
   ax[4].plot(ntime[4:],moving_average(norm(neutron_data[names[i+24]]),5),label=names[i+24],color=colors[i],alpha=0.9)

ax[4].legend()
ax[4].set_ylabel('monitor data shape')

ax[5].plot(xtime,symh)
ax[5].set_ylabel("SymH")

ax[6].plot(xtime,AeInd)
ax[6].set_ylabel("Ae Index")

ax[7].plot(xtime,BZ)
ax[7].set_ylabel("BZ")

plt.xticks(rotation=45)
plt.savefig("dip_vs_neutrons.png")
plt.close()
