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
    return np.convolve(x, np.ones(w), 'same') / w

def norm(x):
    return (x-x.min())/(x.max()-x.min())

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
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()

    return 

def get_pressure_data(file):

    times = []
    field  = []

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

def plot_pressure(time, x, field, ofile, width=5, height=5, nlines=10):

    fig, ax = plt.subplots(figsize=(width,height))

    cmap = plt.cm.cividis(np.linspace(0,1,nlines))

    ## Plot of neutrons
    for i in range(nlines):
        ax.semilogy(time,field[:,i*len(x)//nlines],color=cmap[i],label=f"x:{x[i*len(x)//nlines]}")
    ax.set_ylabel("$\\int \\int |\\boldsymbol{B}|^2 {\\rm d}y {\\rm d}z$")
    ax.set_xlabel("time")

    box = ax.get_position()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#    fig.suptitle(f"Pressure -- {time[0].strftime('%b')} ${time[0].strftime('%Y')}$",size="large")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()

    return 

def get_neutron_data(filename):

    with open(filename, "r") as datFile:
        times = []
        nmdb_data = []
        firstline=0
        for line in datFile:
            # skip all header lines
            if not line.startswith("#"):
                # get station names as first non-commend str
                if firstline==0:
                    stations = line.split()
                    firstline = 1
                else:
                    line_data = []
                    nowrite = 0
                    for i,entry in enumerate(line.split(";")):
                        if entry == "\n":
                            nowrite = 1
                        elif i==0: # record the time
                            times.append(entry)
                        else: # record entry as data
                            if entry.strip(" ").strip("\n") == "null":
                                line_data.append(float("nan"))
                            else:
                                line_data.append(float(entry.strip(" ")))
                    if nowrite == 0:
                        nmdb_data.append(line_data)
    times = pd.to_datetime(times,format="%Y-%m-%d %H:%M:%S")
    return (stations, times, nmdb_data)

def plot_neutrons(stations, times, neutrons, ofile, width=5, height=5):

    fig, ax = plt.subplots(figsize=(width,height))

    cmap = plt.colormaps['cividis']
    colors = cmap(np.linspace(0, 1, len(stations)))

    ## Plot of neutrons
    for i,station in enumerate(stations):
        # Slicing to remove boundaries
        ax.plot(times,moving_average(norm(neutrons[:,i]),6),alpha=0.3,color=colors[i])
        if i == 0:
          mean = np.nan_to_num(moving_average(norm(neutrons[:,i]),6))
        else:
          mean += np.nan_to_num(moving_average(norm(neutrons[:,i]),6))

    # Slicing to remove boundaries
    ax.plot(times,mean/len(stations),color="black",alpha=1)
    ax.set_ylabel("Relative Neutron Counts")
    ax.set_xlim(times.min(),times.max())

#    fig.suptitle(f"Neutron ground data -- {times[0].strftime('%b')} ${times[0].strftime('%Y')}$",size="large")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()

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
    plt.savefig(ofile)
    plt.close()

    return 


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
    plt.savefig(ofile)
    plt.close()

    return

def plot_omni_ae(data,ofile,width=5,height=5):

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

    ax[1,1].plot(data["time"],data["aeind"],color="black")
    ax[1,1].set_ylabel("$\\rm{Ae}$-index")
    ax[1,1].text(0.05, 0.85, "(d)", transform=ax[1,1].transAxes)#,size=20)

    ax[1,0].tick_params('x',labelrotation=45)
    ax[1,1].tick_params('x',labelrotation=45)

#    fig.suptitle(f"Omni data -- {data['time'][0].strftime('%b')} ${data['time'][0].strftime('%Y')}$")
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()

    return



if __name__ == "__main__":
# For nice picture sizes
  pt = 1./72.27 # 72.27 points to an inch.
  jour_sizes = {"APJ": {"onecol": 246.*pt, "twocol": 510.*pt},
                "misc": {"onecol": 246.*pt, "twocol": 372.*pt},
                # Add more journals below or just edit the above numbers
               }

  fig_width = jour_sizes["APJ"]["twocol"]
# Our figure's aspect ratio
  golden = (1 + 5 ** 0.5) / 2
# In figsize - (my_width, my_width/golden)


  for name in ["Aug2018", "May2024", "Mar2022"]:
    (times, x, data) = get_flux_data(f"../../data/{name}/TA16/flux.lst")
    plot_flux(times, x, data, f"../{name}_flux.png",width=fig_width,height=fig_width/golden)
    #(times, x, data) = get_pressure_data(f"../../data/{name}/TA16/pressure.lst")
    #plot_pressure(times, x, data, f"../{name}_pressure.png",width=fig_width,height=fig_width/golden)
    (stations, times, data) = get_neutron_data(f"../../data/{name}/nmdb/nmdb_data_{name}.lst")
    plot_neutrons(stations, times, np.array(data), f"../{name}_neutron_monitor_data.png",width=fig_width,height=fig_width/golden)
    data = get_x_point_data(f"../../data/{name}/TA16/x-point_location_for_all_rates.lst")
#   man_data = get_x_point_data(f"../../data/{name}/TA16/manual_x-point_location.lst")
    plot_x_point(data["time"],data["total"],data["dist"],f"../{name}_x-point_location.png",width=fig_width,height=fig_width/golden)
    data = get_omni_data(f"../../data/{name}/omni/{name}_TA16_parameters.lst")
    plot_omni(data, f"../{name}_omni_data.png",width=fig_width,height=fig_width/golden)
    plot_omni_ae(data, f"../{name}_omni_data_ae.png",width=fig_width,height=fig_width/golden)
