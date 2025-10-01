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

def get_flux_data(file):

    times = []
    field  = []
    x = []

    for i,line in enumerate(np.genfromtxt(file,dtype=None,skip_header=3)):
        times.append(str(line[0]) + '/' + str(line[1]) + '/' + str(line[2]) + '/' + str(line[3]))
        arr = []
        for j in range(4,len(line)):
            if line[j] == "**********":
                arr.append(float('nan'))
            else:
                arr.append(float(line[j]))
        field.append(arr)

        time = pd.to_datetime(times,format="%Y/%j/%H/%M")
        x = [-12 + i for i in range(len(arr))]

    return (time, np.array(x), np.array(field))

def plot_flux(time, x, field, ofile, width=5, height=5, nlines=10):

    fig, ax = plt.subplots(figsize=(width,height))

    cmap = plt.cm.cividis(np.linspace(0,1,nlines))

    ## Plot of neutrons
    handles = []
    labels = []
    for i in range(nlines):
        ax.plot(time,field[:,i*len(x)//nlines],color=cmap[i])[0]
        handles.append(matplotlib.lines.Line2D([],[],color=cmap[i],marker=".",linestyle=""))
        labels.append(f"${x[i*len(x)//nlines]}$")
    ax.set_ylabel("$\\int \\int |B_{\\rm x}| {\\rm d}y {\\rm d}z$")
    ax.set_xlabel("Time")

    box = ax.get_position()
    leg = ax.legend(handles,labels,
                    loc='upper left', bbox_to_anchor=(1, 1.05),title="$\\underline{x}$")
    leg_title = leg.get_title()
    leg_title.set_position( (32,0) )
#    for l in leg.legend_handles:
#      l.set_linestyle("")
#      l.set_marker(".")

#    fig.suptitle(f"Flux -- {time[0].strftime('%b')} ${time[0].strftime('%Y')}$",size="large")
    plt.xticks(rotation=45)
    plt.show()

    return 

if __name__=="__main__":

  for name in ["Mar2022"]:#, "Feb2022", "Jun2015", "May2024", "Oct2024", "Mar2022"]:
    (times,x,data) = get_flux_data(f"../data/{name}/TA16/flux.lst")
    plot_flux(times,x,data,"asdfg.png")
