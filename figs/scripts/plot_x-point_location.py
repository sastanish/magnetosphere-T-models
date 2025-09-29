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

    return {"time":time, "rate":rate, "dist":dist}

if __name__ == "__main__":

  pt = 1./72.27 
  jour_sizes = {"APJ": {"onecol": 246.*pt, "twocol": 510.*pt},
                "misc": {"onecol": 246.*pt, "twocol": 372.*pt},
               }

  width = jour_sizes["APJ"]["twocol"]
  golden = (1 + 5 ** 0.5) / 2

  fig, ax = plt.subplots(nrows=3,ncols=2,figsize=(2*width*0.8,0.8*3*width/golden),dpi=100)
  cmap = plt.cm.cividis(np.linspace(0,1,10))

  ax[0,0].set_title('$|\\boldsymbol{\\Sigma^*}|$')
  ax[0,1].set_title('Reconnection Location')

  for i,name in enumerate(["Aug2018", "Mar2022", "May2024"]):
    TS05_data = get_x_point_data(f"../../data/{name}/TS05/x-point_location.lst")
    TA16_data = get_x_point_data(f"../../data/{name}/TA16/x-point_location.lst")

    ax[i,0].plot(TA16_data["time"],TA16_data["rate"],linestyle="",marker=".",label="rate",alpha=0.3,color=cmap[2])
    ax[i,0].plot(TA16_data["time"],moving_average(TA16_data["rate"],6),linestyle="-",marker="",label="avg rate",color=cmap[2])
    ax[i,0].plot(TS05_data["time"],TS05_data["rate"],linestyle="",marker="x",label="rate",alpha=0.3,color=cmap[8])
    ax[i,0].plot(TS05_data["time"],moving_average(TS05_data["rate"],6),linestyle="--",marker="",label="avg rate",color=cmap[8])
    ax[i,0].set_yscale('log')

    ax[i,0].set_ylabel(TA16_data["time"][0].strftime("%B %Y"))

    ax[i,1].plot(TA16_data["time"],TA16_data["dist"],linestyle="",marker=".",alpha=0.5,color=cmap[2])
    ax[i,1].plot(TS05_data["time"],TS05_data["dist"],linestyle="",marker="x",alpha=0.5,color=cmap[8])
    ax[i,0].tick_params(axis="x",rotation=45)
    ax[i,1].tick_params(axis="x",rotation=45)

    ax[i,1].set_ylim( (1,12) )


  ax[0,0].text(0.02, 0.9, "(a)", transform=ax[0,0].transAxes)#,size=20)
  ax[0,1].text(0.02, 0.9, "(b)", transform=ax[0,1].transAxes)#,size=20)
  ax[1,0].text(0.02, 0.9, "(c)", transform=ax[1,0].transAxes)#,size=20)
  ax[1,1].text(0.02, 0.9, "(d)", transform=ax[1,1].transAxes)#,size=20)
  ax[2,0].text(0.02, 0.9, "(e)", transform=ax[2,0].transAxes)#,size=20)
  ax[2,1].text(0.02, 0.9, "(f)", transform=ax[2,1].transAxes)#,size=20)

  leg = ax[0,1].legend(
     [matplotlib.lines.Line2D([],[],color=cmap[2],marker=".",linestyle="-" ),
      matplotlib.lines.Line2D([],[],color=cmap[8],marker="x",linestyle="--")],
     ["$\\mathcal{M}_1$","$\\mathcal{M}_2$"],
     loc="lower right"
     )

  plt.tight_layout()
  plt.savefig(f"../all_x-point_locations_comparison.png")
  plt.close()
