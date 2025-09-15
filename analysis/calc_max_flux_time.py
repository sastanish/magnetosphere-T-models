import numpy as np
import xarray as xr
import pandas as pd
import datetime

def get_flux_data(file):

    times = []
    field  = []
    x = []

    for i,line in enumerate(np.genfromtxt(file,dtype=None,skip_header=3)):
        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        arr = []
        for j in range(4,len(line)):
            if line[j] == "**********":
                arr.append(float('nan'))
            else:
                arr.append(float(line[j]))
        field.append(arr)

        time = pd.to_datetime(times,format="%Y/%j/%H/%M")
        x = [-12 + i for i in range(len(arr))]

    return (time, np.array(x), np.array(field))

if __name__ == "__main__":

  for name in ["Aug2018", "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    (times, x, data) = get_flux_data(f"../data/{name}/TA16/flux.lst")
    print("---------------" + name + "----------------")
    for ind, xi in enumerate(x):
      max_ind = np.argmax(data[:,ind])
      print("Max flux at x=" + str(xi) + " is " + str(data[max_ind,ind]) + " at " + str(times[max_ind]) )
