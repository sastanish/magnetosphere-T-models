import numpy as np
import pandas as pd
import xarray as xr
import TsyganenkoWrapper as TS
import time as td
from multiprocessing import Pool

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
    stime = td.time()

    [bx, by, bz] = TS.compute.run_ta16(parmod,ps,x,y,z)
    [in_bx, in_by, in_bz] = TS.compute.run_igrf_dipole(iyear,iday,ihour,imin,isec,vgsex,vgsey,vgsez,x,y,z)

    bx += in_bx
    by += in_by
    bz += in_bz

    rate = np.zeros( (nx, ny, nz) )
    rate = TS.compute.total_rate(x,y,z,bx,by,bz)

    time = str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M")).replace(" ","_")

    # Create and write dataset
    ds = xr.Dataset(data_vars={
                          "bx": (["x", "y", "z"], bx),
                          "by": (["x", "y", "z"], by),
                          "bz": (["x", "y", "z"], bz),
                          "rate": (["x", "y", "z"], rate[:,:,:]),
                          },
                   coords={
                              "x": x,
                              "y": y,
                              "z": z,
                      },
                    attrs={"time":time}
               )
    ds.to_netcdf(f"multi_file_{time}.nc",mode='w',format="NETCDF4", 
                 engine='h5netcdf', encoding={
                     "bx":{"zlib":True, "complevel": 7},
                     "by":{"zlib":True, "complevel": 7},
                     "bz":{"zlib":True, "complevel": 7},
                     "rate":{"zlib":True, "complevel": 7}}
                     )
    print(f"some proc just finished, elapsed: {td.time() - stime}")
    return

if __name__ == '__main__':

    # Read data
    omni_data = np.genfromtxt('input_data.lst',dtype=None)

    # Setup the desired GSW Coordinates and data-structure
    (nx, ny, nz) = (50, 50, 40)
    x0 = -10
    x1 = 0
    y0 = -5
    y1 = 5
    z0 = -4
    z1 = 4

    x = np.linspace(x0, x1, nx)
    y = np.linspace(y0, y1, ny)
    z = np.linspace(z0, z1, nz)

    # Set up process pool and operate the compute function on each entry
    # in the omni_data array, or a subset of the array.
    StartTime = td.time()
    with Pool(10) as comp_pool:
        comp_pool.map(compute,omni_data)
    EndTime = td.time()

    print(f"Time Elapsed: {EndTime-StartTime}")
