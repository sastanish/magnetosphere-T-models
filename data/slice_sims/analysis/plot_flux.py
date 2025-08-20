import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import datetime
import scipy

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

    return (time, x, np.array(field))

def plot(time, x, field, ofile, width=5, nlines=10):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    cmap = plt.cm.Blues(np.linspace(0,1,nlines))

#    exp = lambda t,a,b,c,d: a*np.exp(b*(t-d))+c
#    ((a,b,c,d), misc) = scipy.optimize.curve_fit(exp,  x,  np.average(field,axis=0),  p0=(1, 1, -3, -100 ), maxfev = 100000)

    ## Plot of neutrons
    for i in range(nlines):
        ax.plot(time,field[:,i*len(x)//nlines],color=cmap[i],label=f"x:{x[i*len(x)//nlines]}")
    ax.set_ylabel("$\\int \\int \\vec{B}_x dx dy$")
    ax.set_xlabel("time")
    ax.legend()

    fig.suptitle("Flux for storm: " + times[0].strftime("%Y - %m"),size="large")
    plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

  matplotlib.use("AGG")

  for sid in [f"s0{i+1}" for i in range(8)]:
    (times, x, data) = get_field_data(f"../{sid}/flux_Bx.lst")
    plot(times, x, data, f"../figs/{sid}/flux.png")
