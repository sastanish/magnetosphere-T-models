import xarray as xr
import numpy as np
import pandas as pd
import multiprocessing
import re
import os

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

def compute(path,niters=7):

    #load
    ds = xr.open_dataset(path).sel(x=slice(-15,-1))

    #calc pressure and rate
    pressure = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)
    rate = ds.rate

    (cradi,crate) = compute_via_critical_points(pressure.sel(y=0,method="nearest"),rate.sel(y=0,method="nearest"))
    #loop through num iters
    closestPoint = 50
    for j in range(niters):
        try:
            ball = find_ball(pressure,'min',rad=1)
            rateBall = rate.where(pressure == ball)

            loc = rateBall.argmax(dim=("x","y","z"))
            maxrate = rateBall.isel(loc)
            r = np.sqrt(maxrate.x**2+maxrate.y**2+maxrate.z**2)

            if r < closestPoint:
                closestPoint = r
                best_rate = maxrate

            #Artificially increase region to find new one
            pressure += (ball + 1e6)
        except:
            break

    time = str(ds.attrs["time"])

    return (time, float(best_rate), float(closestPoint), crate, cradi)

if __name__ == '__main__':

    dates = ["2018-08-25", "2022-03-13", "2023-03-22", "2024-03-03", "2024-08-11", "2021-11-03", "2022-10-22", "2023-04-23", "2024-03-24", "2024-10-10", "2022-01-14", "2023-02-26", "2023-11-06", "2024-05-10"]
    Nproc = 10

    for date in dates:
        directory = "/users/xnb22215/magnetosphere-T-models/data/TA16/" + date

        filenames = []
        pattern = re.compile(r"output_data_.*\.nc")
        for fname in os.listdir(directory):
            if pattern.match(fname):
                filenames.append(str(directory + "/" + fname))

        with multiprocessing.Pool(Nproc) as pool:
            output = pool.map(compute,filenames)

        header = '''Format:
           1) Time
           2) Maximal reconnection rate at closest pressure dip
           3) Distance from earth
           4) Rate from critical
           5) Distance from critical\n'''
        with open( directory + "/x-point_locations.txt", "w") as f:
            f.write(header)
            for line in output:
                f.write( line[0] + '   ' )
                f.write( f"{line[1]:.4f}" + '   ' )
                f.write( f"{line[2]:.4f}" + '   ' )
                f.write( f"{line[3]:.4f}" + '   ' )
                f.write( f"{line[4]:.4f}" + '\n' )
