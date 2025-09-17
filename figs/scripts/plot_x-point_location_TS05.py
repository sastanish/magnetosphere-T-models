import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import datetime
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size': 16})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
matplotlib.use('AGG')

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
    ax[0].set_ylabel('$|\\boldsymbol{\\Sigma}|$')
    ax[0].set_yscale('log')

    ax[1].plot(time,dist,color="tab:red",linestyle="",marker=".")
    ax[1].set_ylabel('reconnection location')

    ax[1].set_ylim( (1,12) )

    plt.xticks(rotation=45)
    plt.tight_layout()
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

  for name in ["Aug2018", "May2024", "Jun2015", "May2024", "Oct2024"]:
    data = get_x_point_data(f"../../data/{name}/TS05/x-point_location.lst")
    plot_x_point(data["time"],data["rate"],data["dist"],f"../{name}_x-point_location_TS05.png")
