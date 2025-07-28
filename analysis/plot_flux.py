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

    for line in np.genfromtxt(file,dtype=None):
        print(type(line))
        print(len(line))

        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        arr = []
        for i in range(4,len(line)):
            if line[i] == "********":
                arr.append(float('nan'))
            else:
                arr.append(float(line[i]))
        field.append(arr)

    time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return (time, np.array(field))

def plot(time, field, ofile, width=4):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    ## Plot of neutrons
    ax.imshow(field)
    ax.set_ylabel("date")
    ax.set_xlabel("x")

    fig.suptitle("Flux for storm: " + times[0].strftime("%Y - %m"),size="large")
    #plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

    (times, data) = get_field_data("../data/july_scans/s06/flux_Bx.lst")
    plot(times, data, "../figs/july/s06/flux.png")
