import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.style as mplstyle
import multiprocessing

def find_min_ball(ds,ballSize=0.5):

    minima = ds.min()
    loc  = ds.where(ds==ds.min(),drop=True)
    loc = (float(loc.x),float(loc.y),float(loc.z))
    ball = ds.where( np.sqrt( (ds.x - loc[0])**2 + (ds.y - loc[1])**2 + (ds.z - loc[2])**2) <= ballSize, drop=True )
    
    return ball

def find_max_ball(ds,ballSize=0.1):

    maxima = ds.max()
    loc  = ds.where(ds==ds.max(),drop=True)
    loc = (float(loc.x),float(loc.y),float(loc.z))
    ball = ds.where( np.sqrt( (ds.x - loc[0])**2 + (ds.y - loc[1])**2 + (ds.z - loc[2])**2) <= ballSize, drop=True )
    
    return ball

def critical_points(ds,cutoff=10):
    ds = ds.where(ds>cutoff)
    signChanges = np.ones(ds.shape)
    for c in ds.coords:
        if ds[c].size > 1:
            diff = ds.differentiate(c)
            signChanges *= np.sign(diff.shift({c:1})) == np.sign(diff.shift({c:-1}))
            signChanges *= np.sign(diff.shift({c:1})) != np.sign(diff)

    return ds.where( signChanges )

def diagnostic_plot(ds,ball,i):

    mplstyle.use('fast')
    fig, ax = plt.subplots(figsize = (12,6))

    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    pressure.plot.imshow(x="x",y="z",ax=ax,alpha=0.2)

    rate.where( pressure == ball ).plot.imshow(x="x",y="z",ax=ax,cmap="inferno")
    ds.plot.streamplot(x="x",y="z",u="bx",v="bz",ax=ax,color="white")

    canvas = FigureCanvas(fig)
    canvas.print_figure(f"./figs/{i:04d}.png")
    plt.close(fig)

    return 

def compute(ind,niters=3):

    #load
    ds = xr.open_dataset(f"./2024_data_{ind+1}.nc").sel(x=slice(-8,-1))

    #calc pressure and rate
    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    #loop through num iters
    closestPoint = 50
    for j in range(niters):
        #find ball around min pressure
        ball = find_min_ball(pressure)
        rateBall = rate.where(pressure == ball)

        loc = rateBall.argmax(dim=("x","y","z"))
        (xp,yp,zp) = (ds.isel(loc).x,ds.isel(loc).y,ds.isel(loc).z)

        r = float(np.sqrt(xp**2 + yp**2 + zp**2))
        
        if r < closestPoint:
            closestPoint = r
            best_rate = rateBall.max()
            avg_rate = find_max_ball(rateBall).mean()
            xpoint = [float(xp),float(yp),float(zp)]
            
        #Remove last point
        pressure = pressure.where(pressure != ball)

    return (float(best_rate), float(avg_rate), xpoint[0], xpoint[1], xpoint[2])

def compute2(ind):

    #load
    ds = xr.open_dataset(f"../data/2024_data_{ind+1}.nc").sel(x=slice(-8,-1))

    #calc pressure and rate
    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    ball = find_min_ball(pressure)
    rateBall = rate.where(pressure == ball)

    rate = rateBall.max()
    loc = rateBall.argmax(dim=("x","y","z"))
    (xp,yp,zp) = (ds.isel(loc).x,ds.isel(loc).y,ds.isel(loc).z)
    r = float(np.sqrt(xp**2 + yp**2 + zp**2))

    return (float(rate), r)



if __name__ == '__main__':

    Nproc = 20
    names = [1+i for i in range(1199)]

    with multiprocessing.Pool(Nproc) as pool:
        output = pool.map(compute2,names)
    np.savetxt('./x-point_rate.txt',output,fmt="%.4f",header=' Format: \n   1) Maximal reconnection rate at closest pressure dip \n  2) Average around maximal reconnection rate \n  3) x-position \n   4) y-position \n   5) z-position')
