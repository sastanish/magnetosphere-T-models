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

def get_omni_data(file):

  SymHc  = []
  Nind = []
  BZ = []
  avg_BZ = []
  times = []

  for line in np.genfromtxt(file,dtype=None):

    times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

  time = pd.to_datetime(times)

  return {"time":time}


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

  Storms = ["Aug2018", "May2024", "Feb2022", "Jun2015", "Oct2024", "Mar2022"]
  Itimes = [(222,280), (17,199),  (300,313), ( 43, 91), ( 55,187), (  1,103)]
  for name,inds in zip(Storms,Itimes):
    for ind in inds:
      rate  = xr.open_dataarray(f"../../data/{name}/TS05/" + f"rate_{ind}.nc").isel(y=8)
      field = xr.open_dataset(f"../../data/{name}/TS05/" + f"output_{ind}.nc").isel(y=8)

      rate.load()
      field.load()

      rate = rate.where(field.x**2+field.z**2 >= 0.9)
      field = field.where(field.x**2+field.z**2 >= 0.9)

      omni_data = get_omni_data(f"../../data/{name}/omni/{name}_TS05_parameters.lst")

      #pressure = field.bx**2 + field.by**2 + field.bz**2

      #levels = np.geomspace(np.min(pressure),np.max(pressure)/10,num=10)
      #rate_max = rate.where(rate==rate.max(),drop=True)
      #x_center = float(rate_max.x)
      #z_center = float(rate_max.z)
      #print(x_center)
      #points = [ (x_center - 0.5 + i*0.1, z_center - 0.5 + i*0.1) for i in range(10)]

      fig, ax = plt.subplots(figsize = (fig_width,fig_width/golden))
 
      rate.plot.imshow(x="x",y="z",ax=ax,cmap="inferno",vmax=10)
      #pressure.plot.contour(x="x",y="z",levels=levels,ax=ax,cmap="cividis")
      #field.plot.streamplot(x="x",y="z",u="bx",v="bz",ax=ax,color="tab:red",linewidth=0.2,broken_streamlines=False,arrowstyle="<-",arrowsize=0.2,density=1,start_points=points)
      field.plot.streamplot(x="x",y="z",u="bx",v="bz",ax=ax,color="tab:red",linewidth=0.2,broken_streamlines=False,arrowstyle="<-",arrowsize=0.2,density=1)

      earth = plt.Circle((0,0), 1, color="tab:blue", zorder=2.01)
      ax.add_patch(earth)

      ax.set_xlim((-12,-0.01))
      ax.set_ylim((-6,6))
      ax.set_xlabel("x")
      ax.set_ylabel("z")
      plt.tight_layout()
      ax.set_title("")
      #plt.savefig(f"../{name}_slices/{name}_slice_{ind}.png")
      time_str = omni_data["time"][ind].strftime("%Y-%m-%d_%H-%M")
      plt.savefig(f"../{name}_slice_TS05_{time_str}.png")
      plt.close()
      rate.close()
      field.close()
