import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd


def get_symh(file):

  symh = []
  time = []

  for line in np.genfromtxt(file,dtype=None):

    symh.append(float(line[23]))
    time.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

  return {"time":pd.to_datetime(time,format="%Y/%j/%H/%M"), "symh":np.array(symh)}

def get_xpoint(file):

  time = []
  rate = []
  dist = []

  for line in np.genfromtxt(file,dtype=None,skip_header=6):

    time.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
    rate.append( float(line[4]))
    dist.append( float(line[5]))

  return {"time":pd.to_datetime(time,format="%Y/%j/%H/%M"), "rate":rate, "dist":dist}

def get_data_at_symh_min():

  symh_min = []
  symh_min_time = []
  rate_at_min = []
  dist_at_min = []

  for data_dir in [f"../s0{i+1}/" for i in range(8)]:

    omni_data = get_symh(data_dir + 'input_data.lst')
    point_data = get_xpoint(data_dir + 'x-point_location.lst')

    ind = np.argmin(omni_data["symh"])
    min_date = omni_data["time"][ind]

    c = 0
    while min_date != point_data['time'][c]:
      c += 1

    symh_min.append(omni_data['symh'][c])
    symh_min_time.append(min_date)
    rate_at_min.append(point_data['rate'][c])
    dist_at_min.append(point_data['dist'][c])

  return {
          "time":symh_min_time, 
          "symh":symh_min,
          "rate":rate_at_min,
          "dist":dist_at_min
          }


if __name__ == "__main__":

  matplotlib.use("AGG")

  data = get_data_at_symh_min()

  width = 3
  fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(2*width*1.61803,width))

  ax[0].plot(data["symh"],data["dist"],linestyle="",marker="x")
  ax[0].set_ylabel("x-point distance")
  ax[0].set_xlabel("SymH")
  ax[1].plot(data["symh"],data["rate"],linestyle="",marker="x",label=data["time"])
  ax[1].set_ylabel("Reconnection Rate")
  ax[1].set_xlabel("SymH")

  plt.tight_layout()
  plt.savefig("../figs/min_symh_vs_min_x-point.png")
  plt.close()
