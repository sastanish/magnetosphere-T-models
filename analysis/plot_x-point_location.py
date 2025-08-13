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

def get_x_point_data(file):

    times = []
    rate  = []
    dist  = []

    for line in np.genfromtxt(file,dtype=None,skip_header=6):

        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        rate.append( float(line[4]))
        dist.append( float(line[5]))

    time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return {"time":time, "rate":rate, "dist":dist}

def plot_x_point(time,rate,dist,ofile,width=4):

    fig, ax = plt.subplots(nrows=2,figsize=(width*1.61803,2*width),sharex=True)

    ax[0].plot(time,rate,color="tab:blue",linestyle="",marker=".",label="rate")
    ax[0].set_ylabel('rate')
    ax[0].set_yscale('log')

    ax[1].plot(time,dist,color="tab:red",linestyle="",marker=".")
    ax[1].set_ylabel('x-point distance')

    fig.suptitle("x-point location for storm: " + time[0].strftime("%Y - %m"),size="large")
    plt.xticks(rotation=45)
    plt.savefig(ofile)
    plt.close()

    return 

if __name__ == "__main__":

  for sid in ["s01", "s02", "s06"]:
    try:
      data = get_x_point_data(f"../data/{sid}/x-point_location.lst")
      plot_x_point(data["time"],data["rate"],data["dist"],f"../figs/{sid}/x-point_location.png")
    except:
      print("Issue with Storm ID: " + sid)
