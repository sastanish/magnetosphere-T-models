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

    for line in np.genfromtxt(file,dtype=None,skip_header=6):

        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        if line[4] == "********":
            field.append(float('nan'))
        else:
            field.append(float(line[4]))

    time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return (time, field)

def plot(time, field, ofile, width=4):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    ## Plot of neutrons
    ax.plot(time,field)
    ax.set_ylabel("|B|")

    fig.suptitle("Tailward average field for storm: " + times[0].strftime("%Y - %m"),size="large")
    plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

    for sid in ["s01", "s02", "s06"]:
      try:
        (times, data) = get_field_data(f"../data/{sid}/field_strength.lst")
        plot(times, data, f"../figs/{sid}/field_strength.png")
      except:
        print("Issue with Storm ID: " + sid)
