import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import datetime

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_field_data(file):

    times = []
    field  = []
    x = []

    for i,line in enumerate(np.genfromtxt(file,dtype=None)):
        if i==0:
            x.append(float(line[4][2:]))
            for j in range(5,len(line)):
                x.append(float(line[j]))
        else: 
            times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
            arr = []
            for j in range(4,len(line)):
                if line[j] == "**********":
                    arr.append(float('nan'))
                else:
                    arr.append(float(line[j]))
            field.append(arr)

        time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return (time, np.array(x), np.array(field))

def plot(time, x, field, ofile, width=5, nlines=20):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    cmap = plt.cm.Blues(np.linspace(0,1,nlines))

    ## Plot of neutrons
    for i in range(nlines):
        ax.semilogy(time,field[:,i*len(x)//nlines],color=cmap[i],label=f"x:{x[i*len(x)//nlines]}")
    ax.set_ylabel("$\\int \\int |\\vec{B}|^2 dx dy$")
    ax.set_xlabel("time")
    ax.legend()

    fig.suptitle("pressure for storm: " + times[0].strftime("%Y - %m"),size="large")
    plt.xticks(rotation=45)
    fig.tight_layout()
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

  matplotlib.use("AGG")

  for sid in [f"s0{i+1}" for i in range(8)]:
    (times, x, data) = get_field_data(f"../{sid}/pressure.lst")
    plot(times, x, data, f"../figs/{sid}/pressure.png")
