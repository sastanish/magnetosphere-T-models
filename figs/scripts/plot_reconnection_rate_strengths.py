import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import os
import glob

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size': 16})
#rc('text', usetex=True)
#rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
matplotlib.use('AGG')

def get_times(omni_file,num_lines,start_line=0):

  times = []
  lines = np.genfromtxt(omni_file,dtype=None)

  for line in lines[start_line:start_line+num_lines]:
    times.append( str(line[0])+"/"+str(line[1])+"/"+str(line[2])+"/"+str(line[3]) )

  return pd.to_datetime(times,format="%Y/%j/%H/%M")

root_dir = os.environ["PROJ_ROOT"]

pt = 1./72.27 
jour_sizes = {"APJ": {"onecol": 246.*pt, "twocol": 510.*pt},
              "misc": {"onecol": 246.*pt, "twocol": 372.*pt}
             }
fig_width = jour_sizes["APJ"]["twocol"]
golden = (1 + 5 ** 0.5) / 2

for date in ["May2024"]:#,"Mar2022"]:

  fig, ax = plt.subplots(nrows=2,figsize=(fig_width,2*fig_width/golden))

  for i,model in enumerate(["TA16","TS05"]):
    print(i)

    directory = root_dir + "/data/" + date + "/" + model

    print(directory)
    ds = xr.open_mfdataset(directory + "/rates_*.nc",concat_dim=["t"],combine="nested")
    ds = ds.isel(y=8).where(ds.x<=-2).where(ds.x>=-11.5).where(np.abs(ds.z)<5.5)

    times = get_times(directory + "/input_data.lst",len(ds.t))

    ax[i].plot(times,ds["Lorrenz rate"].max(dim=["x","z"],skipna=True),linestyle="",marker=".",label="Lorrenz")
    ax[i].plot(times,ds["Winding rate"].max(dim=["x","z"],skipna=True),linestyle="",marker=".",label="Winding")
#    ax[i].plot(times,ds["Aligned rate"].max(dim=["x","z"],skipna=True),linestyle="",marker=".",label="Aligned")

    ax[i].set_ylim(0.1,500)

  ax[0].legend()

  plt.xticks(rotation=45)
  plt.tight_layout()
  plt.savefig(date+"_rates.png")
