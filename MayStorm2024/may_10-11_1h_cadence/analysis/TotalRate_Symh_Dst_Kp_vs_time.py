import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import xarray as xr
import pandas as pd

def calculate_total_rate(filename):
    rates = [0,0]
    for m,model in enumerate(['TA16', 'TS05']):
        # total reconnection rate
        for i,term in enumerate(["c2_t1", "c2_t2", "c2_t3"]):
            ds = xr.open_dataset(f"../{model}/metrics/" + term + "_" + filename ,engine = "netcdf4",decode_cf=False,decode_times=False)
            if i == 0:
                rate = ds[term].where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 1).sel(y=0,method='nearest')
            else:
                rate += ds[term].where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 1).sel(y=0,method='nearest')
            ds.close()
        rates[m] = rate.max()
    return rates
 
if __name__ == '__main__':
    symh = []
    kp = []
    dst = []
    ap = []
    total_rate = []
    times = []

    # Read in data sources
    TS05_input_data = np.genfromtxt("./rates_data/TS05_data_file.dat",dtype=None)
    omni_index_data = np.genfromtxt("./rates_data/omni_1h_avg_indicies.dat")

    rate_files = []
    with open('filenames.txt') as nfile:
        for line in nfile.readlines():
            rate_files.append(line[0:-1])

    # Fill arrays with data
    for i in range(len(rate_files)):
        times.append(pd.to_datetime(str(TS05_input_data[i][0]) + "_" + str(TS05_input_data[i][1]) + "_" + str(TS05_input_data[i][2]) + "_" + str(TS05_input_data[i][3]), format="%Y_%j_%H_%M"))
        symh.append(TS05_input_data[i][12])
        kp.append(omni_index_data[i,3])
        dst.append(omni_index_data[i,4])
        ap.append(omni_index_data[i,5])
        total_rate.append(calculate_total_rate(rate_files[i]))

    total_rate = np.array(total_rate)
    time_ax = pd.DatetimeIndex(times) # Convert to full pandas time index for nice plots

    # Now we graph
    mplstyle.use('fast')
    fig, axs = plt.subplots(nrows=6,ncols=1,figsize=(8,6*4),sharex=True)

    axs[0].semilogy(time_ax,total_rate[:,0])
    axs[0].set_ylabel("TA16 Max Rate")
    axs[1].semilogy(time_ax,total_rate[:,1])
    axs[1].set_ylabel("TS05 Max Rate")
    axs[2].plot(time_ax,symh)
    axs[2].set_ylabel("SymH")
    axs[3].plot(time_ax,kp)
    axs[3].set_ylabel("Kp Index")
    axs[4].plot(time_ax,ap)
    axs[4].set_ylabel("Ap Index")
    axs[5].plot(time_ax,dst)
    axs[5].set_ylabel("DST Index")

    canvas = FigureCanvas(fig)
    canvas.print_figure("./TotalRate_Symh_Dst_Kp_vs_time.png")
    plt.close(fig)


