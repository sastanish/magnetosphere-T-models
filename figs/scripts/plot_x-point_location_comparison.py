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

#  name = "May2024"
#  name = "Aug2018"
  name = "Mar2022"
  TS05_data = get_x_point_data(f"../../data/{name}/TS05/x-point_location.lst")
  TA16_data = get_x_point_data(f"../../data/{name}/TA16/x-point_location.lst")

  fig, ax = plt.subplots(nrows=2,figsize=(fig_width,2*fig_width/golden),sharex=True)

  cmap = plt.cm.cividis(np.linspace(0,1,10))

  ax[0].plot(TA16_data["time"],TA16_data["rate"],linestyle="",marker=".",label="rate",alpha=0.3,color=cmap[2])
  ax[0].plot(TA16_data["time"],moving_average(TA16_data["rate"],6),linestyle="-",marker="",label="avg rate",color=cmap[2])
  ax[0].plot(TS05_data["time"],TS05_data["rate"],linestyle="",marker="x",label="rate",alpha=0.3,color=cmap[8])
  ax[0].plot(TS05_data["time"],moving_average(TS05_data["rate"],6),linestyle="--",marker="",label="avg rate",color=cmap[8])
  ax[0].set_ylabel('$|\\boldsymbol{\\Sigma^*}|$')
  ax[0].set_yscale('log')
  ax[0].text(0.02, 0.92, "(a)", transform=ax[0].transAxes)#,size=20)

  ax[1].plot(TA16_data["time"],TA16_data["dist"],linestyle="",marker=".",alpha=0.5,color=cmap[2])
  ax[1].plot(TS05_data["time"],TS05_data["dist"],linestyle="",marker="x",alpha=0.5,color=cmap[8])
  ax[1].set_ylabel('Reconnection Location')
  ax[1].text(0.02, 0.92, "(b)", transform=ax[1].transAxes)#,size=20)

  ax[1].set_ylim( (1,12) )

  leg = ax[1].legend(
     [matplotlib.lines.Line2D([],[],color=cmap[2],marker=".",linestyle="-" ),
      matplotlib.lines.Line2D([],[],color=cmap[8],marker="x",linestyle="--")],
     ["TA16","TS05"],
     loc="lower right"
     )

  plt.xticks(rotation=45)
  plt.tight_layout()
  plt.savefig(f"../{name}_x-point_location_comparison.png")
  plt.close()
