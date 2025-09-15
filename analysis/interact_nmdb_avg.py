import numpy as np
import pandas as pd
import datetime
import matplotlib
matplotlib.use('qtagg')
import matplotlib.pyplot as plt

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w



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

def plot_neutrons(stations, times, neutrons, ofile, width=5, height=5):

    fig, ax = plt.subplots(figsize=(width,height))

    cmap = plt.colormaps['cividis']
    colors = cmap(np.linspace(0, 1, len(stations)))

    ## Plot of neutrons
    for i,station in enumerate(stations):
        # Slicing to remove boundaries
        ax.plot(times,moving_average(norm(neutrons[:,i]),10),alpha=0.3,color=colors[i])
        if i == 0:
          mean = np.nan_to_num(norm(neutrons[:,i]))
        else:
          mean += np.nan_to_num(norm(neutrons[:,i]))

    # Slicing to remove boundaries
    ax.plot(times,mean/len(stations),color="black",alpha=1)
    ax.set_ylabel("Relative Neutron Counts")
    ax.set_xlim(times.min(),times.max())

#    fig.suptitle(f"Neutron ground data -- {times[0].strftime('%b')} ${times[0].strftime('%Y')}$",size="large")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    return 


if __name__=="__main__":

  for name in ["Aug2018"]:#, "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    (stations, times, data) = get_neutron_data(f"../data/{name}/nmdb/nmdb_data_{name}.lst")
    plot_neutrons(stations, times, np.array(data), f"../{name}_neutron_monitor_data.png")
