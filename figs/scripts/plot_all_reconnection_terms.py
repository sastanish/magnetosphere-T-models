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
    total_rate  = []
    Lorr_rate  = []
    Wind_rate  = []
    Alig_rate  = []
    dist  = []

    for line in np.genfromtxt(file,dtype=None,skip_header=6):

        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        total_rate.append( float(line[4]))
        Lorr_rate.append( float(line[5]))
        Wind_rate.append( float(line[6]))
        Alig_rate.append( float(line[7]))
        dist.append( float(line[8]))

    time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return {"time":time, "total":total_rate, "lorr":Lorr_rate, "wind":Wind_rate, "alig":Alig_rate, "dist":dist}

if __name__ == "__main__":

  pt = 1./72.27 
  jour_sizes = {"APJ": {"onecol": 246.*pt, "twocol": 510.*pt},
                "misc": {"onecol": 246.*pt, "twocol": 372.*pt},
               }

  width = jour_sizes["APJ"]["twocol"]
  golden = (1 + 5 ** 0.5) / 2

  fig, ax = plt.subplots(nrows=3,ncols=2,figsize=(2*width*0.8,0.8*3*width/golden),dpi=100)
  cmap = plt.cm.cividis(np.linspace(0,1,10))

  cmap = plt.cm.cividis(np.linspace(0,1,3))
  vmap = plt.cm.viridis(np.linspace(0,1,4))
  colors = ["black", cmap[1], vmap[1], cmap[2]]

  ax[0,0].set_title('$\\mathcal{M}_1$')
  ax[0,1].set_title('$\\mathcal{M}_2$')

  for i,name in enumerate(["Aug2018", "Mar2022", "May2024"]):
    TS05_data = get_x_point_data(f"../../data/{name}/TS05/x-point_location_for_all_rates.lst")
    TA16_data = get_x_point_data(f"../../data/{name}/TA16/x-point_location_for_all_rates.lst")

    ax[i,0].plot(TA16_data["time"],moving_average(TA16_data["total"],6),linestyle="-",marker="",color=colors[0])
    ax[i,0].plot(TA16_data["time"],moving_average(TA16_data["lorr"] ,6),linestyle="-",marker="",color=colors[1])
    ax[i,0].plot(TA16_data["time"],moving_average(TA16_data["wind"] ,6),linestyle="-",marker="",color=colors[2])
    ax[i,0].plot(TA16_data["time"],moving_average(TA16_data["alig"] ,6),linestyle="-",marker="",color=colors[3])

    ax[i,0].set_yscale('log')
    ax[i,1].set_yscale('log')
    ax[i,0].set_ylabel(TA16_data["time"][0].strftime("%B %Y"))

    ax[i,1].plot(TS05_data["time"],moving_average(TS05_data["total"],6),linestyle="-",marker="",color=colors[0])
    ax[i,1].plot(TS05_data["time"],moving_average(TS05_data["lorr"] ,6),linestyle="-",marker="",color=colors[1])
    ax[i,1].plot(TS05_data["time"],moving_average(TS05_data["wind"] ,6),linestyle="-",marker="",color=colors[2])
    ax[i,1].plot(TS05_data["time"],moving_average(TS05_data["alig"] ,6),linestyle="-",marker="",color=colors[3])

    ax[i,0].tick_params(axis="x",rotation=45)
    ax[i,1].tick_params(axis="x",rotation=45)



  leg = ax[1,1].legend(
     [matplotlib.lines.Line2D([],[],color=colors[i],linestyle="-") for i in range(4)],
     ["$|\\boldsymbol{\\Sigma^{*}}|$","$|\\boldsymbol{\\Sigma^{*}_{1}}|$","$|\\boldsymbol{\\Sigma^{*}_{2}}|$","$|\\boldsymbol{\\Sigma^{*}_{3}}|$"]
     )

  ax[0,0].text(0.02, 0.9, "(a)", transform=ax[0,0].transAxes)#,size=20)
  ax[0,1].text(0.02, 0.9, "(b)", transform=ax[0,1].transAxes)#,size=20)
  ax[1,0].text(0.02, 0.9, "(c)", transform=ax[1,0].transAxes)#,size=20)
  ax[1,1].text(0.02, 0.9, "(d)", transform=ax[1,1].transAxes)#,size=20)
  ax[2,0].text(0.02, 0.9, "(e)", transform=ax[2,0].transAxes)#,size=20)
  ax[2,1].text(0.02, 0.9, "(f)", transform=ax[2,1].transAxes)#,size=20)

  leg = ax[1,1].legend(
     [matplotlib.lines.Line2D([],[],color=colors[i],linestyle="-") for i in range(4)],
     ["$|\\boldsymbol{\\Sigma^{*}}|$","$|\\boldsymbol{\\Sigma^{*}_{1}}|$","$|\\boldsymbol{\\Sigma^{*}_{2}}|$","$|\\boldsymbol{\\Sigma^{*}_{3}}|$"]
     )

  plt.tight_layout()
  plt.savefig(f"../all_reconnection_rates_comparison.png")
  plt.close()
