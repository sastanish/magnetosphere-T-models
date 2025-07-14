import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_my_omni_data(file):

    SymH   = []
    SymHc  = []
    Nind = []
    avg_BX = []
    avg_BY = []
    avg_BZ = []
    BX = []
    BY = []
    BZ = []
    VX = []
    VY = []
    VZ = []
    tilt = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymH.append(float(line[14]))
        SymHc.append(float(line[23]))
        Nind.append(float(line[22]))
        avg_BX.append(float(line[19]))
        avg_BY.append(float(line[20]))
        avg_BZ.append(float(line[21]))
        BX.append(float(line[4]))
        BY.append(float(line[5]))
        BZ.append(float(line[6]))
        VX.append(float(line[7]))
        VY.append(float(line[8]))
        VZ.append(float(line[9]))
        tilt.append(float(line[17]))

    time = pd.to_datetime(times)

    return [SymH, SymHc, Nind, avg_BX, avg_BY, avg_BZ, BX, BY, BZ, VX, VY, VZ, tilt, times]


def get_tsyganenko_omni_data(file):

    SymH   = []
    SymHc  = []
    Nind = []
    avg_BX = []
    avg_BY = []
    avg_BZ = []
    VX = []
    VY = []
    VZ = []
    tilt = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymH.append(float(line[12]))
        SymHc.append(float(line[19]))
        Nind.append(float(line[17]))
        avg_BX.append(float(line[4]))
        avg_BY.append(float(line[5]))
        avg_BZ.append(float(line[6]))
        VX.append(float(line[7]))
        VY.append(float(line[8]))
        VZ.append(float(line[9]))
        tilt.append(float(line[15]))

    time = pd.to_datetime(times)

    return  [SymH, SymHc, Nind, avg_BX, avg_BY, avg_BZ, VX, VY, VZ, tilt, times]

if __name__ == "__main__":

    mydata = get_my_omni_data("./data/omni_5min_2018_08_25_TA16_parameters.lst")
    tsdata = get_tsyganenko_omni_data("../../data/tsyganenko_data/2018_OMNI_5m_with_RBF_TA16_drivers.dat")

    fig, ax = plt.subplots(nrows=10,figsize=(8,40),sharex=True)

    ax[0].plot(mydata[-1],mydata[0],color="tab:red")
    (xm,xp) = ax[0].get_xlim()
    ax[0].plot(tsdata[-1],tsdata[0],color="tab:blue")
    ax[0].set_ylabel('Symh')

    ax[1].plot(mydata[-1],mydata[1],color="tab:red")
    ax[1].plot(tsdata[-1],tsdata[1],color="tab:blue")
    ax[1].set_ylabel('Symhc')

    ax[2].plot(mydata[-1],mydata[2],color="tab:red")
    ax[2].plot(tsdata[-1],tsdata[2],color="tab:blue")
    ax[2].set_ylabel('N Index')

    ax[3].plot(mydata[-1],mydata[3],color="tab:red")
    ax[3].plot(mydata[-1],mydata[6],color="tab:red",linestyle="",marker=".",alpha=0.3)
    ax[3].plot(tsdata[-1],tsdata[3],color="tab:blue")
    ax[3].set_ylabel('Bx')

    ax[4].plot(mydata[-1],mydata[4],color="tab:red")
    ax[4].plot(mydata[-1],mydata[7],color="tab:red",linestyle="",marker=".",alpha=0.3)
    ax[4].plot(tsdata[-1],tsdata[4],color="tab:blue")
    ax[4].set_ylabel('By')

    ax[5].plot(mydata[-1],mydata[5],color="tab:red")
    ax[5].plot(mydata[-1],mydata[8],color="tab:red",linestyle="",marker=".",alpha=0.3)
    ax[5].plot(tsdata[-1],tsdata[5],color="tab:blue")
    ax[5].set_ylabel('Bz')

    ax[6].plot(mydata[-1],mydata[9],color="tab:red")
    ax[6].plot(tsdata[-1],tsdata[6],color="tab:blue")
    ax[6].set_ylabel('Vx')

    ax[7].plot(mydata[-1],mydata[10],color="tab:red")
    ax[7].plot(tsdata[-1],tsdata[7],color="tab:blue")
    ax[7].set_ylabel('Vy')

    ax[8].plot(mydata[-1],mydata[11],color="tab:red")
    ax[8].plot(tsdata[-1],tsdata[8],color="tab:blue")
    ax[8].set_ylabel('Vz')

    ax[9].plot(mydata[-1],mydata[12],color="tab:red")
    ax[9].plot(tsdata[-1],tsdata[9],color="tab:blue")
    ax[9].set_ylabel('Tilt')

    ax[0].set_xlim((xm,xp))

    plt.savefig("./my_inds_vs_Tsyganenko_2018.png")
    plt.close()
