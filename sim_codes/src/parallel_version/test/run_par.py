import numpy as np
import pandas as pd
import xarray as xr
import TsyganenkoWrapper as TS
import time as td

def compute(line):
    ## This function takes in a line of data from the input file and computes
    ## the resulting field via the fortran wrapper.
    parmod = [0 for i in range(10)]
    parmod[0] = line[12] # Pdyn
    parmod[1] = line[23] # <SymHc>
    parmod[2] = line[22] # N-index
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

    tstart = td.time()
    [bx, by, bz] = TS.compute.run_ta16(parmod,ps,x,y,z)
    tend = td.time()
    print(f'elapsed: {tend-tstart}')

    # Create and write dataset
    ds = xr.Dataset(data_vars={
                          "bx": (["x", "y", "z"], bx),
                          "by": (["x", "y", "z"], by),
                          "bz": (["x", "y", "z"], bz),
                          },
                   coords={
                              "x": x,
                              "y": y,
                              "z": z,
                      }
               )
    ds.to_netcdf(f"./ds.nc",mode='w',format="NETCDF4", 
                 engine='h5netcdf', encoding={
                     "bx":{"zlib":True, "complevel": 7},
                     "by":{"zlib":True, "complevel": 7},
                     "bz":{"zlib":True, "complevel": 7}}
                     )

    return

if __name__ == '__main__':

    # Read data
    omni_data = np.genfromtxt('input_data.lst',dtype=None)

    # Setup the desired GSW Coordinates and data-structure
    (nx, ny, nz) = (75, 100, 20)
    x0 = -15
    x1 = 0
    y0 = -10
    y1 = 10
    z0 = -2
    z1 = 2

    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    z = np.linspace(z0, z1, nz)

    # Set up process pool and operate the compute function on each entry
    # in the omni_data array, or a subset of the array.
    compute(omni_data[0])
