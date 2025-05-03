import numpy as np
import pandas as pd
import xarray as xr
import TsyganenkoWrapper as TS
from multiprocessing import Pool

## NOTE ##
# We use the multiprocessing library so that this code can be run in parallel. 
# All of the computation is done in serial, parallization only allows is to 
# compute many different times at once

Nproc = 2

# The fortran codes are loaded through the TsyganenkoWrapper module.
##########

# Read data

omni_data = np.genfromtxt('selected_may_storm_times.dat',dtype=None)

# Setup the desired GSW Coordinates and data-structure

(nx, ny, nz) = (10*20, 8*20, 8*20)
x0 = -10
x1 = 0
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
    parmod = [0 for i in range(10)]
    parmod[0] = line[-4] # Pdyn
    parmod[1] = line[-1] # <SymHc>
    parmod[2] = line[-2] # N-index
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
    modelNumber = 16 # TA16 model
    dipoleNumber = 2 # IGRF_GSW_08 dipole

    ## DOCS ##
    # TS has only one module, compute. This contains two routines, field and metrics.
    # field: computes the total field of the earth at a given time. If modelNumber=5, use TS05 for
    # external field. If modelNumber=16, use TA16.
    # The dipole model of the earth is either dipoleNumber=1 => DIP_08 or dipoleNumber=2 => IGRF
    # Make sure to pass in the correctly formatted parmod for your chosen model (see TA16_RBF.f or TS04c.f for details).
    (bx, by, bz) = TS.compute.field(x,y,z, (iyear, iday, ihour, imin, isec), (vgsex, vgsey, vgsez), parmod, ps, modelNumber, dipoleNumber)

    # Reconnection Metrics
    # Pass in a grid and field and this method will return the following 
    # calculated quantities along the first dimension;
    # [jx, jy, jz, fx, fy, fz, bfpx, bfpy, bfpz, alpha, lambda, c2_t1, c2_t2, c2_t3]
    output_metrics = np.zeros( (14, nx, ny, nz) )
    output_metrics = TS.compute.metrics(x,y,z,bx,by,bz)

    time = str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")).replace(" ","_")

    # Create and write dataset
    ds = xr.Dataset(data_vars={
                          "bx": (["x", "y", "z"], bx),
                          "by": (["x", "y", "z"], by),
                          "bz": (["x", "y", "z"], bz),
                          "c2_t1": (["x", "y", "z"], output_metrics[11,:,:,:]),
                          "c2_t2": (["x", "y", "z"], output_metrics[12,:,:,:]),
                          "c2_t3": (["x", "y", "z"], output_metrics[13,:,:,:]),
                          # add other fields from output_metrics here
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
                      # Specify extra fields encoding here.
                     )

    print(time)

    return

# Set up process pool and operate the compute function on each entry
# in the omni_data array, or a subset of the array.
with Pool(Nproc) as comp_pool:
    comp_pool.map(compute,omni_data)
