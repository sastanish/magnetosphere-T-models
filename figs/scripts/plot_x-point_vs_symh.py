from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import datetime

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
    avg = moving_average(rate,6)

    return {"time":time, "rate":rate, "dist":dist, "avg rate":avg}

def get_omni_data(file):

    SymHc  = []
    Nind = []
    BZ = []
    avg_BZ = []
    flow = []
    times = []
    aeind = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[14]))
        Nind.append(float(line[22]))
        aeind.append(float(line[13]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))
        flow.append(np.sqrt(float(line[7])**2 + float(line[8])**2 + float(line[9])**2) )

    time = pd.to_datetime(times)

    return {"time":time, "Symh":SymHc, "Nind":Nind, "aeind":aeind, "BZ":BZ, "avg_BZ":avg_BZ, "flow":flow}

if __name__ == "__main__":
# For nice picture sizes
  pt = 1./72.27 # 72.27 points to an inch.
  jour_sizes = {"APJ": {"onecol": 246.*pt, "twocol": 510.*pt},
                "misc": {"onecol": 246.*pt, "twocol": 372.*pt},
                # Add more journals below or just edit the above numbers
               }

  fig_width = jour_sizes["APJ"]["twocol"]
# Our figure's aspect ratio
  golden = (1 + 5 ** 0.5) / 2
# In figsize - (my_width, my_width/golden)

  dates = ["Jun2015", "Aug2018", "Feb2022", "Mar2022", "May2024", "Oct2024"]
  times = []
  symh = []
  dist = []
  rate = []
  speed = []
  handles = []

  cmap = plt.cm.cividis(np.linspace(0,1,len(dates)))
  for i,name in enumerate(dates):
    omni_data = get_omni_data(f"../../data/{name}/omni/{name}_TA16_parameters.lst")
    point_data = get_x_point_data(f"../../data/{name}/TA16/x-point_location.lst")

    min_ind = np.argmin(point_data["dist"])

    times.append(omni_data["time"][min_ind].strftime("%b %Y, %H:%M"))
    symh.append(omni_data["Symh"][min_ind])
    speed.append(omni_data["flow"][min_ind])
    dist.append(point_data["dist"][min_ind])
    rate.append(point_data["avg rate"][min_ind])

    handles.append(matplotlib.lines.Line2D([],[],color=cmap[i],marker=".",linestyle=""))

  fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(fig_width,2*fig_width/golden),sharex=True)

  ax[0].scatter(dist,symh,c=cmap)
  ax[1].scatter(dist,speed,c=cmap)

  leg = ax[0].legend(handles,times,ncols=2,
                     mode="expand", bbox_to_anchor=(0, 1.05, 1, 0.3),markerscale=3)

  ax[0].set_ylabel("SymH")
  ax[1].set_ylabel("$|\\boldsymbol{V}|$")
  ax[1].set_xlabel("$R_{\\oplus,\\min}$")

  plt.tight_layout()
  plt.savefig("../x-point_vs_symh_and_v.png")
  plt.close()


