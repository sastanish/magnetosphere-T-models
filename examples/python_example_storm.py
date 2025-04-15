import numpy as np
import pandas as pd
import xarray as xr
import Tsyganenko_wrapper as wrap
from multiprocessing import Pool

## NOTE ##
# We use the multiprocessing library so that this code can be run in parallel. 
# All of the computation is done in serial, parallization only allows is to 
# compute many different times at once

Nproc = 4

# The fortran codes are loaded through the Tsyganenko_wrapper module.
##########

# Read data

omni_data = np.genfromtxt('may_2024_storm_with_TS05_vars.dat',dtype=None)

# Setup the desired GSW Coordinates and data-structure

(nx, ny, nz) = (400, 100, 100)
x0 = -20
x1 = 2
y0 = -4
y1 = 4
z0 = -4
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
    ds.to_netcdf(f"data/{time}.nc")#,mode='w',format="NETCDF4", 
#                engine='h5netcdf', encoding={
#                    "bx":{"zlib":True, "complevel": 7},
#                    "by":{"zlib":True, "complevel": 7},
#                    "bz":{"zlib":True, "complevel": 7}}
#                    )

    print(time)

    # If you wish to calculate the reconnection metrics, uncomment
    # the following function call

#    recon_metrics(ds)

    return

compute(omni_data[1])
# Set up process pool and operate the compute function on each entry
# in the omni_data array, or a subset of the array.
#with Pool(Nproc) as comp_pool:
#    comp_pool.map(compute,omni_data[0:2])
