import numpy as np
import xarray as xr
import pandas as pd
import datetime

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

def get_omni_data(file):

    SymHc  = []
    aeind = []
    BZ = []
    avg_BZ = []
    flow = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[14]))
        aeind.append(float(line[13]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))
        flow.append(np.sqrt(float(line[7])**2 + float(line[8])**2 + float(line[9])**2) )

    time = pd.to_datetime(times)

    return {"time":time, "Symh":SymHc, "aeind":aeind, "BZ":BZ, "avg_BZ":avg_BZ, "flow":flow}

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



if __name__ == "__main__":

  times = []
  symh = []
  dist = []
  rate = []
  for i,name in enumerate(["Aug2018", "Jun2015", "May2024", "Mar2022"]):
    omni_data = get_omni_data(f"../data/{name}/omni/{name}_TA16_parameters.lst")
    point_data = get_x_point_data(f"../data/{name}/TA16/x-point_location.lst")

    min_ind = np.argmin(point_data["dist"])

    times.append(omni_data["time"][min_ind].strftime("%b %Y at %H:%M"))
    symh.append(omni_data["Symh"][min_ind])
    dist.append(point_data["dist"][min_ind])
    rate.append(point_data["avg rate"][min_ind])

  print(f"""
  \\begin{{center}}
  \\begin{{tabular}}{{ c | c c c c }}
                               & {times[0]} & {times[1]} & {times[2]} & {times[3]} \\\\
  $R_\\oplus$                  &${dist[0]}$ &${dist[1]}$ &${dist[2]}$ &${dist[3]}$ \\\\
  $|\\boldsymbol{{\\Sigma^*}}|$&${rate[0]:.2f}$ &${rate[1]:.2f}$ &${rate[2]:.2f}$ &${rate[3]:.2f}$ \\\\
  SymH                         &${symh[0]}$ &${symh[1]}$ &${symh[2]}$ &${symh[3]}$
  \\end{{tabular}}
  \\end{{center}}
  """)
