import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import multiprocessing

def plot(ind):

    ds = xr.open_dataset(f"../data/2024-03_data_{ind}.nc")

    mplstyle.use('fast')
    fig, ax = plt.subplots(figsize = (12,6))

    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    pressure.sel(y=0,method='nearest').plot.imshow(x="x",y="z",ax=ax,cmap="bwr")

    rate.where( rate > 10 ).sel(y=0,method='nearest').plot.imshow(x="x",y="z",ax=ax,cmap="inferno")

    ds.sel(y=0,method='nearest').plot.streamplot(x="x",y="z",u="bx",v="bz",ax=ax,color="white")

    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./figs/{ind:04d}.png")
    plt.close(fig)
    ds.close()

    return 

if __name__ == '__main__':

    Nproc = 10
    names = [1+i for i in range(283)]

    with multiprocessing.Pool(Nproc) as pool:
        pool.map(plot,names)
