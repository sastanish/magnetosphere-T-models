import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_neutron_data(directory,hskip=42):

    # Clean up neutron monitor data
    neutron_data = np.genfromtxt(directory + "neutron_data.txt",
                                 skip_header=hskip,dtype=None,names=True,missing_values="null"
                                 )
    times = []
    for i in range(len(neutron_data)):
        times.append(pd.to_datetime(str(neutron_data["DATE"][i]) + "_" + str(neutron_data["HOUR"][i]), format="%Y-%m-%d_%H:%M:%S"))

    names = neutron_data.dtype.names
    ntime = pd.to_datetime(times)

    return (names, ntime, neutron_data)

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
    fig, ax = plt.subplots(nrows=5,figsize=(8,16),sharex=True)

    xpoints = get_x_point_data(directory)
    omnidata = get_omni_data(directory)
#    colors = matplotlib.color_sequences["tab20"]
#    neutrondata = get_neutron_data(directory)

    ax[0].plot(xpoints["time"],xpoints["brate"],color="tab:red",linestyle="",marker=".",label="max rate",alpha=0.9)
    ax[0].plot(xpoints["time"],xpoints["crate"],color="tab:blue",linestyle="",marker=".",label="max rate",alpha=0.9)
    ax[0].set_ylabel('rate')
    ax[0].set_yscale('log')

    ax[1].plot(xpoints["time"],xpoints["brad"],color="tab:red",linestyle="",marker=".",alpha=0.9)
    ax[1].plot(xpoints["time"],xpoints["crad"],color="tab:blue",linestyle="",marker=".",alpha=0.9)
    ax[1].set_ylabel('radius')

    ax[2].plot(omnidata["time"],omnidata["SymHc"])
    ax[2].set_ylabel("<SymHc>")

    ax[3].plot(omnidata["time"],omnidata["Nind"])
    ax[3].set_ylabel("N Index")

    ax[4].plot(omnidata["time"],omnidata["BZ"],color="black")
    ax[4].plot(omnidata["time"],omnidata["avg_BZ"],color="tab:red")
    ax[4].set_ylabel("BZ")

    plt.xticks(rotation=45)
    plt.savefig("../figs/TA16/" + date + "/x-point_location.png")
    plt.close()

    return

if __name__ == "__main__":

    dates = ["2018-08-25", "2022-03-13", "2023-03-22", "2024-03-03", "2024-08-11", "2021-11-03", "2022-10-22", "2023-04-23", "2024-03-24", "2024-10-10", "2022-01-14", "2023-02-26", "2023-11-06", "2024-05-10"]
    for date in dates:
        print(date)
        plot_statistics(date)

