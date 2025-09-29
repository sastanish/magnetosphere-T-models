from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import datetime

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size': 16})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
matplotlib.use('AGG')

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

def get_omni_data(file):

    SymHc  = []
    Nind = []
    BZ = []
    avg_BZ = []
    flow = []
    times = []
    aeind = []

    for line in np.genfromtxt(file,dtype=None):

        times.append(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

        SymHc.append(float(line[14]))
        Nind.append(float(line[22]))
        aeind.append(float(line[13]))
        BZ.append(float(line[6]))
        avg_BZ.append(float(line[21]))
        flow.append(np.sqrt(float(line[7])**2 + float(line[8])**2 + float(line[9])**2) )

    time = pd.to_datetime(times)

    return {"time":time, "Symh":SymHc, "Nind":Nind, "aeind":aeind, "BZ":BZ, "avg_BZ":avg_BZ, "flow":flow}

def plot_omni(data,ofile,width=5,height=5):

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(2*width,2*height),sharex=True,dpi=100)

#    ax[1].plot(data["time"],data["Nind"])
#    ax[1].set_ylabel("N Index")
    (lbx,lby) = (0.05,0.9)

    ax[0,0].plot(data["time"],data["BZ"],color="black")
    ax[0,0].set_ylabel("$B_{\\rm z}$")
    ax[0,0].text(0.05, 0.1, "(a)", transform=ax[0,0].transAxes)#,size=20)

    ax[1,0].plot(data["time"],data["flow"],color="black")
    ax[1,0].set_ylabel("$V$")
    ax[1,0].text(0.05, 0.85, "(c)", transform=ax[1,0].transAxes)#,size=20)

    plt.xticks(rotation=45)

    ax[0,1].plot(data["time"],data["Symh"],color="black")
    ax[0,1].set_ylabel("Sym-H")
    ax[0,1].text(0.05, 0.1, "(b)", transform=ax[0,1].transAxes)#,size=20)

    ax[1,1].plot(data["time"],data["Nind"],color="black")
    ax[1,1].set_ylabel("$N$-index")
    ax[1,1].text(0.05, 0.85, "(d)", transform=ax[1,1].transAxes)#,size=20)

    ax[1,0].tick_params('x',labelrotation=45)
    ax[1,1].tick_params('x',labelrotation=45)

#    fig.suptitle(f"Omni data -- {data['time'][0].strftime('%b')} ${data['time'][0].strftime('%Y')}$")
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()

    return

if __name__ == "__main__":
  pt = 1./72.27 
  jour_sizes = {"APJ": {"onecol": 246.*pt, "twocol": 510.*pt},
                "misc": {"onecol": 246.*pt, "twocol": 372.*pt},
               }

  width = jour_sizes["APJ"]["twocol"]
  golden = (1 + 5 ** 0.5) / 2


  for name in ["Aug2018", "May2024", "Mar2022"]:
    data = get_omni_data(f"../../data/{name}/TA16/input_data.lst")
    plot_omni(data, f"../{name}_omni_data.png",width*0.7,0.7*width/golden)
