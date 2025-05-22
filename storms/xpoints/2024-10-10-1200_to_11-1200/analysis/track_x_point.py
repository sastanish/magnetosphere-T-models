import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import multiprocessing

def find_ball(ds,minOrMax,rad=0.5):

    if minOrMax== "max":
        loc = ds.isel(ds.argmax(dim=("x","y","z")))
    elif minOrMax== "min":
        loc = ds.isel(ds.argmin(dim=("x","y","z")))
    else:
        print("Must pass either, 'min' or 'max' as second arg.")
        return
    return ds.where( np.sqrt( (ds.x - loc.x)**2 + (ds.y - loc.y)**2 + (ds.z - loc.z)**2) <= rad, drop=True)

def diagnostic_plot(ds,ball,i):

    mplstyle.use('fast')
    fig, ax = plt.subplots(figsize = (12,6))

    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    pressure.sel(y=0,method='nearest').plot.imshow(x="x",y="z",ax=ax,alpha=0.2)

    rate.where( pressure == ball ).sel(y=0,method='nearest').plot.imshow(x="x",y="z",ax=ax,cmap="inferno")
    ds.sel(y=0,method='nearest').plot.streamplot(x="x",y="z",u="bx",v="bz",ax=ax,color="white")

    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./figs/{i:04d}.png")
    plt.close(fig)

    return 

def compute(ind,niters=3):

    #load
    ds = xr.open_dataset(f"../data/2024_data_{ind+1}.nc").sel(x=slice(-8,-1))

    #calc pressure and rate
    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    #loop through num iters
    closestPoint = 50
    for j in range(niters):
        ball = find_ball(pressure,'min',rad=0.5)
        rateBall = rate.where(pressure == ball)

        loc = rateBall.argmax(dim=("x","y","z"))
        maxrate = rateBall.isel(loc)
        r = np.sqrt(maxrate.x**2+maxrate.y**2+maxrate.z**2)
        
        if r < closestPoint:
            closestPoint = r
            best_rate = maxrate
            #avg_ball = find_ball(rateBall,'max',rad=0.1)
            #avg_rate = avg_ball.mean()
            
        #Remove last point
        pressure = pressure.where(pressure != ball)

    return (float(best_rate), float(closestPoint))

def compute2(ind):

    #load
    ds = xr.open_dataset(f"./2024_data_{ind+1}.nc").sel(x=slice(-8,-1))

    #calc pressure and rate
    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    ball = find_ball(pressure,'min',rad=0.5)
    rateBall = rate.where(pressure == ball)

    loc = rateBall.argmax(dim=("x","y","z"))
    maxrate = rateBall.isel(loc)
    r = np.sqrt(maxrate.x**2+maxrate.y**2+maxrate.z**2)
        
    #diagnostic_plot(ds,ball,ind)

    return (float(maxrate), float(r))

if __name__ == '__main__':

    Nproc = 20
    names = [1+i for i in range(1199)]

    with multiprocessing.Pool(Nproc) as pool:
        output = pool.map(compute,names)
    np.savetxt('./x-point_rate.txt',output,fmt="%.4f",header=' Format: \n   1) Maximal reconnection rate at closest pressure dip \n  2) Distance from earth')
