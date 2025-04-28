import numpy as np
import pandas as pd
import xarray as xr
import Tsyganenko_wrapper as wrap
from multiprocessing import Pool

Nproc = 32

# Read data from truncated MayStorm2024 file.

omni_data = np.genfromtxt('may10_15-24hrs_5min_cadence.dat',dtype=None)

# Setup the desired GSW Coordinates and data-structure
# with a resolution of 15 cells per earth radi

(nx, ny, nz) = (15*30, 8*30, 6*30)
x0 = -17
x1 = -2
y0 = -4
y1 = 4
z0 = -2
z1 = 4

x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)
z = np.linspace(z0, z1, nz)

# Computation Functions

def compute(line):
    ## This function takes in a line of data from the input file and computes
    ## the resulting field via the fortran wrapper.
    parmod = [line[16], line[12], line[5], line[6], line[17], line[18], line[19], line[20], line[21], line[22]]
    ps = line[15]
    vgsex = line[7]
    vgsey = line[8] + 29.78
    vgsez = line[9]
    iyear = line[0]
    iday = line[1]
    ihour = line[2]
    imin = line[3]
    isec = 0

    # External field
    [ex_bx, ex_by, ex_bz] = wrap.models.run_ts05(parmod,ps,x,y,z)
    # Internal field via igrf
    [in_bx, in_by, in_bz] = wrap.models.run_igrf_dipole(iyear,iday,ihour,imin,isec,vgsex,vgsey,vgsez,x,y,z)

    time = str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")).replace(" ","_")

    # Create and write dataset
    ds = xr.Dataset(data_vars={
                          "bx": (["x", "y", "z"], ex_bx + in_bx),
                          "by": (["x", "y", "z"], ex_by + in_by),
                          "bz": (["x", "y", "z"], ex_bz + in_bz),
                          },
                   coords={
                              "x": x,
                              "y": y,
                              "z": z,
                      },
                    attrs={"time":time}
               )
    ds = ds
    ds.where(np.sqrt(ds.x**2+ds.y**2+ds.z**2)>0.5).to_netcdf(f"data/{time}.nc",mode='w',format="NETCDF4", 
                 engine='h5netcdf', encoding={
                     "bx":{"zlib":True, "complevel": 7},
                     "by":{"zlib":True, "complevel": 7},
                     "bz":{"zlib":True, "complevel": 7}}
                     )

    print(time)

    return

# Set up process pool and operate the compute function on each entry
# in the omni_data
with Pool(Nproc) as comp_pool:
    comp_pool.map(compute,omni_data)
