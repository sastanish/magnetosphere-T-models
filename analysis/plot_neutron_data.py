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

def get_neutron_data(filename):

    with open(filename, "r") as datFile:
        times = []
        nmdb_data = []
        firstline=0
        for line in datFile:
            # skip all header lines
            if not line.startswith("#"):
                # get station names as first non-commend str
                if firstline==0:
                    stations = line.split()
                    firstline = 1
                else:
                    line_data = []
                    nowrite = 0
                    for i,entry in enumerate(line.split(";")):
                        if entry == "\n":
                            nowrite = 1
                        elif i==0: # record the time
                            times.append(entry)
                        else: # record entry as data
                            if entry.strip(" ").strip("\n") == "null":
                                line_data.append(float("nan"))
                            else:
                                line_data.append(float(entry.strip(" ")))
                    if nowrite == 0:
                        nmdb_data.append(line_data)
    times = pd.to_datetime(times,format="%Y-%m-%d %H:%M:%S")
    return (stations, times, nmdb_data)

def plot(stations, times, neutrons, ofile, width=4):

    fig, ax = plt.subplots(figsize=(width*1.61803,width))

    ## Plot of neutrons
    for i,station in enumerate(stations):
        ax.plot(times,moving_average(norm(neutrons[:,i]),10),alpha=0.5)
    ax.set_ylabel("Relative Neutron Counts")

    fig.suptitle("neutron ground data for storm: " + times[0].strftime("%Y - %m"),size="large")
    plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

    for sid in ["s01", "s02", "s06"]:
      try:
        (stations, times, data) = get_neutron_data(f"../data/nmdb/nmdb_data_{sid}.lst")
        plot(stations, times, np.array(data), f"../figs/{sid}/neutron_monitor_data.png")
      except:
        print("Issue with Storm ID: " + sid)

