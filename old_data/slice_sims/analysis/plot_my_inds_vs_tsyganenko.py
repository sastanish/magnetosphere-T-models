import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_my_omni_data(file):

    SymHc  = []
    Nind = []
    avg_BZ = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[23]))
        Nind.append(float(line[22]))
        avg_BZ.append(float(line[21]))

    time = pd.to_datetime(times)

    return {"time":time, "SymHc":SymHc, "Nind":Nind, "avg_BZ":avg_BZ}


def get_tsyganenko_omni_data(file):

    SymHc  = []
    Nind = []
    avg_BZ = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[19]))
        Nind.append(float(line[17]))
        avg_BZ.append(float(line[6]))

    time = pd.to_datetime(times)

    return {"time":time, "SymHc":SymHc, "Nind":Nind, "avg_BZ":avg_BZ}

if __name__ == "__main__":

    mydata = get_my_omni_data("../data/prepared_omni_data/TA16/data/2018-08-25_TA16_parameters.lst")
    tsdata = get_tsyganenko_omni_data("../data/tsyganenko_data/2018_OMNI_5m_with_RBF_TA16_drivers.dat")

    fig, ax = plt.subplots(nrows=3,figsize=(8,8),sharex=True)

    ax[0].plot(mydata["time"],mydata["SymHc"],color="tab:red")
    (xm,xp) = ax[0].get_xlim()
    ax[0].plot(tsdata["time"],tsdata["SymHc"],color="tab:blue")
    ax[0].set_ylabel('SymHc')

    ax[1].plot(mydata["time"],mydata["Nind"],color="tab:red")
    ax[1].plot(tsdata["time"],tsdata["Nind"],color="tab:blue")
    ax[1].set_ylabel('N index')

    ax[2].plot(mydata["time"],mydata["avg_BZ"],color="tab:red")
    ax[2].plot(tsdata["time"],tsdata["avg_BZ"],color="tab:blue")
    ax[2].set_ylabel('$<BZ>$')

    ax[0].set_xlim((xm,xp))
    ax[1].set_xlim((xm,xp))
    ax[2].set_xlim((xm,xp))

    plt.savefig("../figs/" + "my_inds_vs_Tsyganenko_2018.png")
    plt.close()
