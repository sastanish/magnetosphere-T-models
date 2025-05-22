import xarray as xr
import numpy as np
import pandas as pd
import multiprocessing

def find_critical_points(ds,dim):

    diff = ds.differentiate(dim)
    right_sign_changes = ds.where(diff * diff.shift({dim:1}) <= 0)
    left_sign_changes = ds.where(diff * diff.shift({dim:-1}) <= 0)

    return right_sign_changes.notnull() + left_sign_changes.notnull()

def find_ball(ds,minOrMax,rad=0.5):

    if minOrMax== "max":
        loc = ds.isel(ds.argmax(dim=("x","y","z")))
    elif minOrMax== "min":
        loc = ds.isel(ds.argmin(dim=("x","y","z")))
    else:
        print("Must pass either, 'min' or 'max' as second arg.")
        return
    return ds.where( np.sqrt( (ds.x - loc.x)**2 + (ds.y - loc.y)**2 + (ds.z - loc.z)**2) <= rad, drop=True)

def compute_via_critical_points(pressure,rate):

    points = rate.where(find_critical_points(pressure,"x")*find_critical_points(pressure,"z"))
    points = points.where(np.sqrt(points.x**2+points.z**2)>1.5)

    stacked_points = points.stack(s=["x","z"]).notnull()
    stacked_rate = rate.stack(s=["x","z"])

    rates = stacked_rate[stacked_points]
    radii = np.sqrt(rates.x**2+rates.z**2)

    try:
        cradi = float(radii.min())
        crate = float(rates[radii.argmin()])
    except:
        cradi = np.nan
        crate = np.nan

    return (cradi, crate)

def compute(ind,niters=3):

    #load
    ds = xr.open_dataset(f"../data/2024-03_data_{ind+1}.nc").sel(x=slice(-8,-1))

    #calc pressure and rate
    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.mag_c2_t1 + ds.mag_c2_t2 + ds.mag_c2_t3

    (cradi,crate) = compute_via_critical_points(pressure.sel(y=0,method="nearest"),rate.sel(y=0,method="nearest"))

    #loop through num iters
    closestPoint = 50
    for j in range(niters):
        ball = find_ball(pressure,'min',rad=1)
        rateBall = rate.where(pressure == ball)

        loc = rateBall.argmax(dim=("x","y","z"))
        maxrate = rateBall.isel(loc)
        r = np.sqrt(maxrate.x**2+maxrate.y**2+maxrate.z**2)

        if r < closestPoint:
            closestPoint = r
            best_rate = maxrate

        #Remove last point
        pressure = pressure.where(pressure != ball)


    return (float(best_rate), float(closestPoint), crate, cradi)

if __name__ == '__main__':

    Nproc = 10
    names = [1+i for i in range(89)]

    with multiprocessing.Pool(Nproc) as pool:
        output = pool.map(compute,names)
    np.savetxt('./x-point_rate.txt',output,fmt="%.4f",header=' Format: \n   1) Maximal reconnection rate at closest pressure dip \n  2) Distance from earth \n  3) Rate from critical \n  4) Distance from critical')
