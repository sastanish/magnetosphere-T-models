import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle

def get_omni_data(file):

  SymHc  = []
  Nind = []
  BZ = []
  avg_BZ = []
  times = []

  for line in np.genfromtxt(file,dtype=None):

    times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

    SymHc.append(float(line[23]))
    Nind.append(float(line[22]))
    BZ.append(float(line[6]))
    avg_BZ.append(float(line[21]))

  time = pd.to_datetime(times)

  return {"time":time, "SymHc":SymHc, "Nind":Nind, "BZ":BZ, "avg_BZ":avg_BZ}

if __name__=="__main__":
  data_directory = "../data/july_scans/s01/"
  fig_directory = "../figs/july/s01/"
  mplstyle.use('fast')

  for i in range(1,603):
    omni_data = get_omni_data(data_directory + "input_data.lst")
    rate = xr.open_dataarray(data_directory + f"rate_{i}.nc").isel(y=151)
    rate.load()
    field = xr.open_dataset(data_directory + f"output_{i}.nc").isel(y=151)
    field.load()
   
    pressure = field.bx**2 + field.by**2 + field.bz**2
   
    pressure = pressure.where(field.x**2+field.z**2 >= 1)
    rate = rate.where(field.x**2+field.z**2 >= 1)

    levels = np.geomspace(np.min(pressure),np.max(pressure)/10,num=10)
   
    fig, ax = plt.subplots(figsize = (5,5))
   
    rate.plot.imshow(x="x",y="z",ax=ax,cmap="inferno",vmax=10)
    pressure.plot.contour(x="x",y="z",levels=levels,ax=ax,cmap="cividis")
   
    ax.set_title(omni_data["time"][i].strftime("%Y/%m/%d - %H:%M"))
    canvas = FigureCanvas(fig)
    canvas.print_figure(fig_directory + f"tail_slice_{i}.png")
    plt.close(fig)
    rate.close()
    field.close()
