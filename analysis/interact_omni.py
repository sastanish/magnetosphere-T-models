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

def plot_omni(data,width=5,height=5):

    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(2*width,2*height),sharex=True)

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
    plt.show()

    return

if __name__=="__main__":

  for name in ["May2024"]:#, "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    data = get_omni_data(f"../data/{name}/omni/{name}_TA16_parameters.lst")
    plot_omni(data)
