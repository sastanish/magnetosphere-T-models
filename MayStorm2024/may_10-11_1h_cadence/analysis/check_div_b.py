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
        ds = xr.open_dataset(f"../{model}/field/" + filename ,engine = "netcdf4",decode_cf=False,decode_times=False)
        divb = ds.bx.differentiate('x') + ds.by.differentiate('y') + ds.bz.differentiate('z')
        divb.sel(y=0,method="nearest").sel(x=slice(-15,-2)).plot.imshow(x="z",y="x",ax=ax[m])

    fig.suptitle("Div B " + str(ds.attrs["time"]))
    ax[0].set_title("TA16")
    ax[1].set_title("TS05")
    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./divB/divB_{filename[0:-3]}.png")
    plt.close(fig)


if __name__ == '__main__':
    mplstyle.use('fast')

    files = []
    with open('filenames.txt') as nfile:
        for line in nfile.readlines():
            files.append(line[0:-1])

    proc_pool = multiprocessing.Pool(12)
    proc_pool.map(plot, files)
