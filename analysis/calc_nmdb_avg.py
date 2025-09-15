import numpy as np
import pandas as pd
import datetime

def norm(x):
    return (x-x.min())/(x.max()-x.min())


def get_neutron_data(filename):

    with open(filename, "r") as datFile:
        times = []
        nmdb_data = []
        firstline=0
        for line in datFile:
            # skip all header lines
            if not line.startswith("#"):
                # get station names as first non-commend str
                if firstline==0:
                    stations = line.split()
                    firstline = 1
                else:
                    line_data = []
                    nowrite = 0
                    for i,entry in enumerate(line.split(";")):
                        if entry == "\n":
                            nowrite = 1
                        elif i==0: # record the time
                            times.append(entry)
                        else: # record entry as data
                            if entry.strip(" ").strip("\n") == "null":
                                line_data.append(float("nan"))
                            else:
                                line_data.append(float(entry.strip(" ")))
                    if nowrite == 0:
                        nmdb_data.append(line_data)
    times = pd.to_datetime(times,format="%Y-%m-%d %H:%M:%S")
    return (stations, times, nmdb_data)

if __name__=="__main__":

  for name in ["Aug2018", "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    (stations, times, data) = get_neutron_data(f"../data/{name}/nmdb/nmdb_data_{name}.lst")
    data = np.array(data)
    for i,station in enumerate(stations):
      # Slicing to remove boundaries
      if i == 0:
        mean = np.nan_to_num(norm(data[:,i]))
      else:
        mean += np.nan_to_num(norm(data[:,i]))
    mean = mean/len(stations)
    with open(f"../data/{name}/nmdb/mean_nmdb_data_{name}.lst","w") as ofile:
      for (time,val) in zip(times,mean):
        ofile.write(str(time)+"    "+str(val)+"\n")
