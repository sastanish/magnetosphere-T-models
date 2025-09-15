import numpy as np
import xarray as xr
import pandas as pd
import datetime

def get_omni_data(file):

    SymHc  = []
    Nind = []
    BZ = []
    avg_BZ = []
    flow = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[14]))
        Nind.append(float(line[22]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))
        flow.append(np.sqrt(float(line[7])**2 + float(line[8])**2 + float(line[9])**2) )

    time = pd.to_datetime(times)

    return {"time":time, "Symh":SymHc, "Nind":Nind, "BZ":BZ, "avg_BZ":avg_BZ, "flow":flow}


if __name__ == "__main__":

  for name in ["Aug2018", "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    data = get_omni_data(f"../data/{name}/omni/{name}_TA16_parameters.lst")
    print("---------------" + name + "----------------")
    min_ind = np.argmax(data["Nind"])
    print("Max N-index" + " is " + str(data["Nind"][min_ind]) + " at " + str(data["time"][min_ind]) )
