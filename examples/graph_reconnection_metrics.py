import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import multiprocessing

def make_plot(ds,var):
    fig, ax = plt.subplots(figsize = (3,6))
    ds[var].sel(y=0,method='nearest').plot.imshow(x='z',y='x',ax=ax)
    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./figs/{time}_{var}.png")
    plt.close(fig)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
mplstyle.use('fast')

def load_and_write_plots(time)

    ds = xr.open_dataset("./{time}.nc",engine = "netcdf4",decode_cf=False,decode_times=False)
    ds = ds.where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 1)

    make_plot(ds,"mag_j")

    ds.close()

    for var in ["j","t1","t2","t3"]:

        ds = xr.open_dataset("./{time}_{var}.nc",engine = "netcdf4",decode_cf=False,decode_times=False)
        ds = ds.where(np.sqrt(ds.x**2 + ds.y**2 + ds.z**2) > 1)

        make_plot(ds,"mag_{var}")

        ds.close()

if __name__ == '__main__':

    omni_data = np.genfromtxt('may_2024_storm_with_TS05_vars.dat',dtype=None)

    time = lambda line : str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

    times = [time(line) for line in omni_data[0:6]]

    proc_pool = multiprocessing.Pool(3)
    proc_pool.map(load_and_write_plots, times)

