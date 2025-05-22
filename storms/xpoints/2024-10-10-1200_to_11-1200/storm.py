import numpy as np
import pandas as pd
import xarray as xr
import TsyganenkoWrapper as TS
from multiprocessing import Pool

## NOTE ##
# We use the multiprocessing library so that this code can be run in parallel. 
# All of the computation is done in serial, parallization only allows is to 
# compute many different times at once

Nproc = 40

# The fortran codes are loaded through the TsyganenkoWrapper module.
##########

# Read data

omni_data = np.genfromtxt('input_data.lst',dtype=None)

# Setup the desired GSW Coordinates and data-structure

(nx, ny, nz) = (8*25, 4*25, 4*25)
x0 = -8
x1 = 0
y0 = -2
y1 = 2
z0 = -2
z1 = 2

x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)
z = np.linspace(z0, z1, nz)

# Computation Functions

def compute(line):
    ## This function takes in a line of data from the input file and computes
    ## the resulting field via the fortran wrapper.
    parmod = [0 for i in range(10)]
    parmod[0] = line[18] # Pdyn
    parmod[1] = line[24] # <SymHc>
    parmod[2] = line[23] # N-index
    parmod[3] = line[20] # <By IMF>
    ps = line[17]
    vgsex = line[7]
    vgsey = line[8] + 29.78
    vgsez = line[9]
    iyear = line[0]
    iday = line[1]
    ihour = line[2]
    imin = line[3]
    isec = 0

    # External field
    [bx, by, bz] = TS.compute.run_ta16(parmod,ps,x,y,z)
    # Internal field via igrf
    [in_bx, in_by, in_bz] = TS.compute.run_igrf_dipole(iyear,iday,ihour,imin,isec,vgsex,vgsey,vgsez,x,y,z)

    bx += in_bx
    by += in_by
    bz += in_bz

    # Reconnection Metrics
    # Pass in a grid and field and this method will return the following 
    # calculated quantities along the first dimension;
    # [jx, jy, jz, fx, fy, fz, bfpx, bfpy, bfpz, alpha, lambda,
    #  mag_c2_t1, c2_t1_x, c2_t1_y, c2_t1_z,
    #  mag_c2_t2, c2_t2_x, c2_t2_y, c2_t2_z,
    #  mag_c2_t3, c2_t3_x, c2_t3_y, c2_t3_z]
    output_metrics = np.zeros( (22, nx, ny, nz) )

    output_metrics = TS.compute.metrics(x,y,z,bx,by,bz)

    time = str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")).replace(" ","_")

    # Create and write dataset
    ds = xr.Dataset(data_vars={
                          "bx": (["x", "y", "z"], bx),
                          "by": (["x", "y", "z"], by),
                          "bz": (["x", "y", "z"], bz),
                          "mag_c2_t1": (["x", "y", "z"], output_metrics[11,:,:,:]),
                          "mag_c2_t2": (["x", "y", "z"], output_metrics[15,:,:,:]),
                          "mag_c2_t3": (["x", "y", "z"], output_metrics[19,:,:,:]),
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
                     "mag_c2_t1":{"zlib":True, "complevel": 7},
                     "mag_c2_t2":{"zlib":True, "complevel": 7},
                     "mag_c2_t3":{"zlib":True, "complevel": 7}}
                      # Specify extra fields encoding here.
                     )

    print(time)

    return

# Set up process pool and operate the compute function on each entry
# in the omni_data array, or a subset of the array.
with Pool(Nproc) as comp_pool:
    comp_pool.map(compute,omni_data)
