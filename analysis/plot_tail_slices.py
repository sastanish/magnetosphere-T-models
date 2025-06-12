import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import multiprocessing
import re
import os

def plot_tail(inp):

    date = inp[0]
    ind = inp[1]
    ds = xr.open_dataset(f"../data/TA16/{date}/output_data_{ind}.nc")
    # Assume ds is a slice in y
    ds = ds.where(np.sqrt(ds.x**2+ds.z**2) > 1).sel(y=0,method="nearest").squeeze()

    mplstyle.use('fast')
    fig, ax = plt.subplots(figsize = (12,6))

    ds.rate.plot.imshow(x="x",y="z",ax=ax,cmap="inferno")
    ds.plot.streamplot(x="x",y="z",u="bx",v="bz",ax=ax,color="white")

    ax.set_title(str(ds.attrs["time"]))
    canvas = FigureCanvas(fig)
    canvas.print_figure(f"../figs/TA16/{date}/tail_slices/{ind}.png")
    plt.close(fig)

    return 

def compute_for_date(date,Nproc):

    inds = []
    pattern = re.compile(r"output_data_.*\.nc")
    for fname in os.listdir(f"../data/TA16/{date}/"):
        if pattern.match(fname):
            inds.append((date,str(fname).strip("output_data_").strip(".nc")))

    with multiprocessing.Pool(Nproc) as pool:
        output = pool.map(plot_tail,inds)


    return 

if __name__ == "__main__":

    dates = ["2018-08-25", "2022-03-13", "2023-03-22", "2024-03-03", "2024-08-11", "2021-11-03", "2022-10-22", "2023-04-23", "2024-03-24", "2024-10-10", "2022-01-14", "2023-02-26", "2023-11-06", "2024-05-10"]
    Nproc = 10
    for date in dates:
        compute_for_date(date,Nproc)

