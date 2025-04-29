import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import multiprocessing

def plot(filename):

    # First plot magnetic field lines and recon rate.
    fig, ax = plt.subplots(ncols=2, figsize = (6,6))

    for m,model in enumerate(['TA16', 'TS05']):
        # field lines
        bds = xr.open_dataset(f"../{model}/field/" + filename ,engine = "netcdf4",decode_cf=False,decode_times=False)
        bds = bds.where(np.sqrt(bds.x**2 + bds.y**2 + bds.z**2) > 1).sel(y=0,method='nearest')

        # total reconnection rate
        for i,term in enumerate(["c2_t1", "c2_t2", "c2_t3"]):
            ds = xr.open_dataset(f"../{model}/metrics/" + term + "_" + filename ,engine = "netcdf4",decode_cf=False,decode_times=False)
            if i == 0:
                rate = ds[term].where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 1).sel(y=0,method='nearest')
            else:
                rate += ds[term].where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 1).sel(y=0,method='nearest')
            ds.close()
        rate.plot.imshow(x='z',y='x',ax=ax[m])
        bds.plot.streamplot(x='z',y="x",u='bz',v='bx',ax=ax[m])
    
    fig.suptitle(str(bds.attrs["time"]))
    ax[0].set_title("TA16")
    ax[1].set_title("TS05")
    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./comparison/field_and_rate_{filename[0:-3]}.png")
    plt.close(fig)

    # Now compare magnitude and vector angle
    fig, ax = plt.subplots(ncols=2, figsize = (6,6))

    ds16 = xr.open_dataset(f"../TA16/field/" + filename ,engine = "netcdf4",decode_cf=False,decode_times=False).sel(y=0,method='nearest')
    ds05 = xr.open_dataset(f"../TS05/field/" + filename ,engine = "netcdf4",decode_cf=False,decode_times=False).sel(y=0,method='nearest')

    mag05 = np.sqrt(ds05.bx**2 + ds05.by**2 + ds05.bz**2) 
    mag16 = np.sqrt(ds16.bx**2 + ds16.by**2 + ds16.bz**2) 
    magdiff = np.abs(mag05-mag16)

    angle = np.arccos((ds16.bx*ds05.bx + ds16.by*ds05.by + ds16.bz*ds05.bz)/(mag16*mag05))

    angle.plot.imshow(x='z',y='x',ax=ax[0])
    magdiff.plot.imshow(x='z',y='x',ax=ax[1])

    fig.suptitle(str(bds.attrs["time"]))
    ax[0].set_title('angle between')
    ax[1].set_title('mag difference')

    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./comparison/angle_and_mag_{filename[0:-3]}.png")
    plt.close(fig)


if __name__ == '__main__':
    mplstyle.use('fast')

    files = []
    with open('filenames.txt') as nfile:
        for line in nfile.readlines():
            files.append(line[0:-1])

    proc_pool = multiprocessing.Pool(12)
    proc_pool.map(plot, files)
