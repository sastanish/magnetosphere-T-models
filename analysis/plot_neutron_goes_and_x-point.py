import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_goes_particles(filename):

    with open(filename, "r") as datFile:
        times = []
        data_labels = ["P1", "P5", "P10", "P30", "P50", "P100", "E_8", "E2_0", "E4_0", "P60", "P500"]
        goes_data = []
        for i,line in enumerate(datFile): 
            if i != 0: # skip first line
                line_data = []
                nowrite = 0
                # remove bad characters
                for bad in ( "(", ")" ):
                    line = line.replace(bad,"")
                # replace arrow with seperator
                line = line.replace(" ->",",")
                for i,entry in enumerate(line.split(",")):
                    if entry == "\n":
                        nowrite = 1
                    elif i==0: # record the time
                        times.append(entry)
                    else: # record entry as data
                        if entry.strip(" ").strip("\n") == "-100000.0":
                            line_data.append(float("nan"))
                        else:
                            line_data.append(float(entry.strip(" ")))
                if nowrite == 0:
                    goes_data.append(line_data)
    times = pd.to_datetime(times,format="%Y-%m-%d %H:%M:%S")
    return (data_labels, times, goes_data)


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

def get_omni_data(directory):

    SymHc  = []
    Nind = []
    BZ = []
    avg_BZ = []
    times = []

    for line in np.genfromtxt(directory + "input_data.lst",dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[23]))
        Nind.append(float(line[22]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))

    time = pd.to_datetime(times)

    return {"time":time, "SymHc":SymHc, "Nind":Nind, "BZ":BZ, "avg_BZ":avg_BZ}

def get_x_point_data(directory):

    times = []
    brate = []
    brad  = []
    crate = []
    crad  = []

    for line in np.genfromtxt(directory + 'x-point_locations.txt',dtype=None,skip_header=6):

        times.append(str(line[0]).replace("_"," "))
        brate.append(float(line[1]))
        brad.append( float(line[2]))
        crate.append(float(line[3]))
        crad.append( float(line[4]))

    time = pd.to_datetime(times)

    return {"time":time, "brate":brate, "brad":brad, "crate":crate, "crad":crad}

def plot_statistics(date):

    directory = "../data/TA16/" + date + "/"
    nmdb_file = "../data/nmdb_data/nmdb_data_" + date + ".lst"
    goes_file = "../data/goes_data/goes_flux_data_" + date + ".lst"
    fig, ax = plt.subplots(nrows=4,figsize=(10,12),sharex=True)

    ## Use helper functions to get data from data dir
    xpoints = get_x_point_data(directory)
    omnidata = get_omni_data(directory)
    (stations, ntimes, neutrons) = get_neutron_data(nmdb_file)
    (goes_ids, gtimes, goesflux) = get_goes_particles(goes_file)

    neutrons = np.array(neutrons)
    goesflux = np.array(goesflux)

    ## Plot of neutrons
    for i,station in enumerate(stations):
        ax[0].plot(ntimes,moving_average(norm(neutrons[:,i]),10),alpha=0.5)
    ax[0].set_ylabel("Relative Neutron Counts")

    ax[1].plot(xpoints["time"],xpoints["crate"],color="tab:blue",linestyle="",marker=".",label="max rate",alpha=0.9)
    ax[1].plot(xpoints["time"],xpoints["brate"],color="tab:red",linestyle="",marker=".",label="max rate",alpha=0.9)
    ax[1].set_ylabel('rate')
    ax[1].set_yscale('log')

    data_inds = (1,2,3,4,5)
    for ind in data_inds:
        ax[2].plot(gtimes,goesflux[:,ind],label=goes_ids[ind])
        #ax[2].plot(gtimes,moving_average(goesflux[:,ind],10),label=goes_ids[ind])
    ax[2].legend()
    ax[2].set_ylabel("Proton Fluxes")

    ax[3].plot(xpoints["time"],xpoints["brad"],color="tab:red",linestyle="",marker=".",alpha=0.9)
    ax[3].plot(xpoints["time"],xpoints["crad"],color="tab:blue",linestyle="",marker=".",alpha=0.9)
    ax[3].set_ylabel('radius')


    ax[0].axvline(gtimes[np.argmax(goesflux[:,1])],linestyle="--",color="grey")
    ax[1].axvline(gtimes[np.argmax(goesflux[:,1])],linestyle="--",color="grey")
    ax[2].axvline(gtimes[np.argmax(goesflux[:,1])],linestyle="--",color="grey")
    ax[3].axvline(gtimes[np.argmax(goesflux[:,1])],linestyle="--",color="grey")

    fig.suptitle("storm: " + date,size="xx-large")
    plt.xticks(rotation=45)
    plt.savefig("../figs/TA16/" + date + "/protons_neutrons_and_x-points.png")
    plt.close()

    return 

if __name__ == "__main__":

    dates = ["2018-08-25", "2022-03-13", "2023-03-22", "2024-03-03", "2024-08-11", "2021-11-03", "2022-10-22", "2023-04-23", "2024-03-24", "2024-10-10", "2022-01-14", "2023-02-26", "2023-11-06", "2024-05-10"]

    for date in dates:
        print(date)
        plot_statistics(date)
