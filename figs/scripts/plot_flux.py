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

    for i,line in enumerate(np.genfromtxt(file,dtype=None,skip_header=3)):
        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        arr = []
        for j in range(4,len(line)):
            if line[j] == "**********":
                arr.append(float('nan'))
            else:
                arr.append(float(line[j]))
        field.append(arr)

        time = pd.to_datetime(times,format="%Y/%j/%H/%M")
        x = [-12 + i for i in range(len(arr))]

    return (time, np.array(x), np.array(field))

def plot(time, x, field, ofile, width=5, nlines=10):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    cmap = plt.cm.Blues(np.linspace(0,1,nlines))

    ## Plot of neutrons
    for i in range(nlines):
        ax.plot(time,field[:,i*len(x)//nlines],color=cmap[i],label=f"x:{x[i*len(x)//nlines]}")
    ax.set_ylabel("$\\int \\int |\\vec{B}_x| dy dz$")
    ax.set_xlabel("time")
    ax.legend()

    fig.suptitle(f"Flux -- {time[0].strftime('%b')} ${time[0].strftime('%Y')}$",size="large")
    plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":
  matplotlib.use('AGG')
  #matplotlib.use('module://matplotlib-backend-kitty')
  plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Helvetica'
  })

  for name in ["Aug2018", "Feb2022", "Jun2015", "May2024", "Oct2024"]:
    (times, x, data) = get_field_data(f"../../data/{name}/TA16/flux.lst")
    plot(times, x, data, f"../{name}_flux.png")
