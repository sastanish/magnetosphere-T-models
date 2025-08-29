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

def plot(data,ofile,width=4):

    fig, ax = plt.subplots(nrows=3,figsize=(width*1.61803,3*width),sharex=True)

    ax[0].plot(data["time"],data["Symh"])
    ax[0].set_ylabel("SymH")

#    ax[1].plot(data["time"],data["Nind"])
#    ax[1].set_ylabel("N Index")

    ax[1].plot(data["time"],data["BZ"],color="tab:red")
    ax[1].set_ylabel("$|\\mathbf{BZ}|$")

    ax[2].plot(data["time"],data["flow"],color="tab:green")
    ax[2].set_ylabel("$|\\mathbf{V}|$")

    plt.xticks(rotation=45)
    fig.suptitle(f"Omni data -- {data['time'][0].strftime('%b')} ${data['time'][0].strftime('%Y')}$")
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()

    return

if __name__ == "__main__":
  matplotlib.use('AGG')
  #matplotlib.use('module://matplotlib-backend-kitty')
  plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'Helvetica'
  })


  for name in ["Aug2018", "Feb2022", "Jun2015", "May2024", "Oct2024"]:
    data = get_omni_data(f"../../data/{name}/omni/{name}_TA16_parameters.lst")
    plot(data, f"../{name}_omni_data.png")

