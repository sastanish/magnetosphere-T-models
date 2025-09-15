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

def get_x_point_data(file):

    times = []
    rate  = []
    dist  = []

    for line in np.genfromtxt(file,dtype=None,skip_header=6):

        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        rate.append( float(line[4]))
        dist.append( float(line[5]))

    time = pd.to_datetime(times,format="%Y/%j/%H/%M")

    return {"time":time, "rate":rate, "dist":dist}

def plot_x_point(time,srate,sdist,mrate,mdist,ofile,width=5,height=5):

    fig, ax = plt.subplots(nrows=2,figsize=(width,2*height),sharex=True)

#    ax[0].plot(time,srate,color="gray",linestyle="",marker=".",alpha=0.5)
    ax[0].plot(time,srate,color="tab:blue",linestyle="",marker=".",label="rate")
    ax[0].set_ylabel('$|\\boldsymbol{\\Sigma^*}|$')
    ax[0].set_yscale('log')
    ax[0].text(0.02, 0.92, "(a)", transform=ax[0].transAxes)#,size=20)

#    ax[1].plot(time,sdist,color="tab:red",linestyle="",marker=".",alpha=0.5)
    ax[1].plot(time,sdist,color="tab:red",linestyle="",marker=".")
    ax[1].set_ylabel('reconnection location')
    ax[1].text(0.02, 0.92, "(b)", transform=ax[1].transAxes)#,size=20)

    ax[1].set_ylim( (1,12) )

#    fig.suptitle(f"X-point Location -- {time[0].strftime('%b')} ${time[0].strftime('%Y')}$")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    return 


if __name__=="__main__":

  for name in ["Aug2018"]:#, "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    sim_data = get_x_point_data(f"../data/{name}/TA16/x-point_location.lst")
    man_data = get_x_point_data(f"../data/{name}/TA16/manual_x-point_location.lst")
    plot_x_point(sim_data["time"],sim_data["rate"],sim_data["dist"],man_data["rate"],man_data["dist"],f"../{name}_x-point_location.png")
