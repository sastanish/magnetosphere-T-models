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

def plot_x_point(time,rate,dist,ofile,width=5,height=5):

    fig, ax = plt.subplots(nrows=2,figsize=(width,2*height),sharex=True)

    b_cmap = plt.cm.cividis(np.linspace(0,1,10))

    ax[0].plot(time,rate,linestyle="",marker=".",label="rate",alpha=0.3,color=b_cmap[2])
    ax[0].plot(time,moving_average(rate,6),color="black",linestyle="-",label="avg rate")
    ax[0].set_ylabel('$|\\boldsymbol{\\Sigma^*}|$')
    ax[0].set_yscale('log')
    ax[0].text(0.02, 0.92, "(a)", transform=ax[0].transAxes)#,size=20)

    ax[1].plot(time,dist,color=b_cmap[1],linestyle="",marker=".")
    ax[1].set_ylabel('Reconnection Location')
    ax[1].text(0.02, 0.92, "(b)", transform=ax[1].transAxes)#,size=20)

    ax[1].set_ylim( (1,12) )

#    fig.suptitle(f"X-point Location -- {time[0].strftime('%b')} ${time[0].strftime('%Y')}$")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    return 

def get_x_point_data(file):

    times = []
    total_rate  = []
    Lorr_rate  = []
    Wind_rate  = []
    Alig_rate  = []
    dist  = []

    for line in np.genfromtxt(file,dtype=None,skip_header=6):

        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        total_rate.append( float(line[4]))
        Lorr_rate.append( float(line[5]))
        Wind_rate.append( float(line[6]))
        Alig_rate.append( float(line[7]))
        dist.append( float(line[8]))

    time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return {"time":time, "total":total_rate, "lorr":Lorr_rate, "wind":Wind_rate, "alig":Alig_rate, "dist":dist}

if __name__=="__main__":

  for name in ["May2024"]:#, "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    data = get_x_point_data(f"../data/{name}/TA16/x-point_location_for_all_rates.lst")
    plot_x_point(data["time"],data["total"],data["dist"],f"../{name}_x-point_location.png")
