import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import datetime
import scipy

def get_omni_data(file):

    SymHc  = []
    Nind = []
    BZ = []
    avg_BZ = []
    tilt = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[23]))
        Nind.append(float(line[22]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))
        tilt.append(float(line[17]))

    time = pd.to_datetime(times)

    return {"time":time, "SymHc":SymHc, "Nind":Nind, "BZ":BZ, "avg_BZ":avg_BZ, "tilt":tilt}


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

def plot(time, x, field, otime, tilt, ofile, width=5, nlines=10):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    cmap = plt.cm.Blues(np.linspace(0,1,nlines))

#    exp = lambda t,a,b,c,d: a*np.exp(b*(t-d))+c
#    ((a,b,c,d), misc) = scipy.optimize.curve_fit(exp,  x,  np.average(field,axis=0),  p0=(1, 1, -3, -100 ), maxfev = 100000)

    ## Normalize and remove shape of tilt
    tilt = (tilt-np.min(tilt))/(np.max(tilt)-np.min(tilt)) * (np.max(field)-np.min(field)) - np.min(field)

    ## Plot of neutrons
    for i in range(nlines):
        ax.plot(time,field[:,i*len(x)//nlines]+tilt[:-1],color=cmap[i],label=f"x:{x[i*len(x)//nlines]}")
    ax.set_ylabel("$\\int \\int \\vec{B}_x dx dy$")
    ax.set_xlabel("time")
    ax.legend()

    fig.suptitle("Flux for storm without tilt: " + times[0].strftime("%Y - %m"),size="large")
    plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

    (times, x, data) = get_field_data("../data/july_scans/s06/flux_Bx.lst")
    omni_data = get_omni_data("../data/july_scans/s06/input_data.lst")
    plot(times, x, data, omni_data["time"], np.array(omni_data["tilt"]), "../figs/july/s06/flux_without_tilt.png")
