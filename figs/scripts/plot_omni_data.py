import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_omni_data(file):

    SymHc  = []
    Nind = []
    BZ = []
    avg_BZ = []
    tilt = []
    times = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[23]))
        Nind.append(float(line[22]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))
        tilt.append(float(line[17]))

    time = pd.to_datetime(times)

    return {"time":time, "SymHc":SymHc, "Nind":Nind, "BZ":BZ, "avg_BZ":avg_BZ, "tilt":tilt}

def plot(data,ofile,width=4):

    fig, ax = plt.subplots(nrows=2,figsize=(width*1.61803,2*width),sharex=True)

    ax[0].plot(data["time"],data["SymHc"])
    ax[0].set_ylabel("<SymHc>")

#    ax[1].plot(data["time"],data["Nind"])
#    ax[1].set_ylabel("N Index")

    ax[1].plot(data["time"],data["BZ"],color="black")
    ax[1].plot(data["time"],data["avg_BZ"],color="tab:red")
    ax[1].set_ylabel("BZ")

#    ax[3].plot(data["time"],data["tilt"],color="tab:green")
#    ax[3].set_ylabel("tilt")

    plt.xticks(rotation=45)
    fig.suptitle("Omni data for storm: " + data["time"][0].strftime("%Y - %m"),size="large")
    plt.savefig(ofile)
    plt.close()

    return

if __name__ == "__main__":
    for sid in ["s01", "s02", "s06"]:
        try:
            data = get_omni_data(f"../data/{sid}/input_data.lst")
            plot(data,f"../figs/{sid}/omni_data.png")
        except:
            print("Issue with Storm ID: " + sid)
        
