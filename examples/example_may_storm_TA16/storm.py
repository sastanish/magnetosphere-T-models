import numpy as np
import pandas as pd
import xarray as xr
import Tsyganenko_wrapper as wrap
from multiprocessing import Pool

## NOTE ##
# We use the multiprocessing library so that this code can be run in parallel. 
# All of the computation is done in serial, parallization only allows is to 
# compute many different times at once

Nproc = 2

# The fortran codes are loaded through the Tsyganenko_wrapper module.
##########

# Read data

omni_data = np.genfromtxt('selected_may_storm_times.dat',dtype=None)

# Setup the desired GSW Coordinates and data-structure

(nx, ny, nz) = (15*10, 8*10, 6*10)
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
    parmod = [0 for i in range(10)]
    parmod[0] = line[-4] # Pdyn
    parmod[1] = line[-1] # <SymHc>
    parmod[2] = line[-2] # I am using the N-index here. Discuss with David.
    parmod[3] = line[5] # <By IMF>
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
    [bx, by, bz] = wrap.models.run_ta16(parmod,ps,x,y,z)
    # Internal field via igrf
    [in_bx, in_by, in_bz] = wrap.models.run_igrf_dipole(iyear,iday,ihour,imin,isec,vgsex,vgsey,vgsez,x,y,z)
    # Full field
    bx = bx+in_bx
    by = by+in_by
    bz = bz+in_bz

    print("trying metrics")

    # Reconnection Metrics
    [c2_t1, c2_t2, c2_t3] = wrap.models.calculate_metrics(x,y,z,bx,by,bz)

    time = str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")).replace(" ","_")

    print("now writing metrics")

    # Create and write dataset
    ds = xr.Dataset(data_vars={
                          "bx": (["x", "y", "z"], bx),
                          "by": (["x", "y", "z"], by),
                          "bz": (["x", "y", "z"], bz),
                          "c2_t1": (["x", "y", "z"], c2_t1),
                          "c2_t2": (["x", "y", "z"], c2_t2),
                          "c2_t3": (["x", "y", "z"], c2_t3),
                          },
                   coords={
                              "x": x,
                              "y": y,
                              "z": z,
                      },
                    attrs={"time":time}
               )
    ds.to_netcdf(f"data/{time}.nc",mode='w',format="NETCDF4", 
                 engine='h5netcdf', encoding={
                     "bx":{"zlib":True, "complevel": 7},
                     "by":{"zlib":True, "complevel": 7},
                     "bz":{"zlib":True, "complevel": 7},
                     "c2_t1":{"zlib":True, "complevel": 7},
                     "c2_t2":{"zlib":True, "complevel": 7},
                     "c2_t3":{"zlib":True, "complevel": 7}}
                     )

    print(time)

    return

# Set up process pool and operate the compute function on each entry
# in the omni_data array, or a subset of the array.
with Pool(Nproc) as comp_pool:
    comp_pool.map(compute,omni_data)
