import numpy as np
import pandas as pd
import xarray as xr
import Tsyganenko_wrapper as wrap
from multiprocessing import Pool

## NOTE ##
# We use the multiprocessing library so that this code can be run in parallel. 
# All of the computation is done in serial, parallization only allows is to 
# compute many different times at once

Nproc = 1

# The fortran codes are loaded through the Tsyganenko_wrapper module.
##########

# Read data

omni_data = np.genfromtxt('may_2024_storm_with_TS05_vars.dat',dtype=None)

# Setup the desired GSW Coordinates and data-structure

(nx, ny, nz) = (256, 64, 64)
x0 = -25
x1 = 2
y0 = -4
y1 = 4
z0 = -4
z1 = 4

x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)
z = np.linspace(z0, z1, nz)

# Computation Functions

def recon_metrics(ds):
## This function computes the reconnection metrics and writes them out.
## Turn off and on inside main computation loop
    mag_B = np.sqrt(ds.bx**2 + ds.by**2 + ds.bz**2)

    jx = ds.bz.differentiate("y") - ds.by.differentiate("z")
    jy = ds.bx.differentiate("z") - ds.bz.differentiate("x")
    jz = ds.by.differentiate("x") - ds.bx.differentiate("y")

    mag_J = np.sqrt(jx**2 + jy**2 + jz**2)

    alpha = (jx * ds.bx + jy * ds.by + jz * ds.bz)/mag_B**2

 
    Fx = jy * ds.bz - jz * ds.by
    Fy = jz * ds.bx - jx * ds.bz
    Fz = jx * ds.by - jy * ds.bx
 
    mag_F = np.sqrt(Fx**2 + Fy**2 + Fz**2)
    
 
    B_fperpx = (ds.by * Fz - ds.bz * Fy)/mag_F 
    B_fperpy = (ds.bz * Fx - ds.bx * Fz)/mag_F 
    B_fperpz = (ds.bx * Fy - ds.by * Fx)/mag_F 

 
    lamb = (jx * B_fperpx + jy * B_fperpy + jz * B_fperpz)/mag_B**2
 
    omega1 = (B_fperpz.differentiate("y") - B_fperpy.differentiate("z")) * Fx/mag_F \
           + (B_fperpx.differentiate("z") - B_fperpz.differentiate("x")) * Fy/mag_F \
           + (B_fperpy.differentiate("x") - B_fperpx.differentiate("y")) * Fz/mag_F
 
    omega2 = (B_fperpz.differentiate("y") - B_fperpy.differentiate("z")) * B_fperpx/mag_B**2\
           + (B_fperpx.differentiate("z") - B_fperpz.differentiate("x")) * B_fperpy/mag_B**2\
           + (B_fperpy.differentiate("x") - B_fperpx.differentiate("y")) * B_fperpz/mag_B**2
 
    C2_t1_coef = (lamb*omega1 - (lamb.differentiate("x")*ds.bx + lamb.differentiate("y")*ds.by + lamb.differentiate("z")*ds.bz) )
    C2_t2_coef = lamb*(alpha + omega2)
    grad_alpha_x = alpha.differentiate("x")
    grad_alpha_y = alpha.differentiate("y")
    grad_alpha_z = alpha.differentiate("z")

 
    C2_t1_x = C2_t1_coef * Fx/mag_F
    C2_t1_y = C2_t1_coef * Fy/mag_F
    C2_t1_z = C2_t1_coef * Fz/mag_F
    mag_C2_t1 = np.sqrt(C2_t1_x**2 + C2_t1_y**2 + C2_t1_z**2)

 
    C2_t2_x = C2_t2_coef * B_fperpx
    C2_t2_y = C2_t2_coef * B_fperpy
    C2_t2_z = C2_t2_coef * B_fperpz
    mag_C2_t2 = np.sqrt(C2_t2_x**2 + C2_t2_y**2 + C2_t2_z**2)

 
    C2_t3_x = grad_alpha_y * ds.bz - grad_alpha_z * ds.by
    C2_t3_y = grad_alpha_z * ds.bx - grad_alpha_x * ds.bz
    C2_t3_z = grad_alpha_x * ds.by - grad_alpha_y * ds.bx
    mag_C2_t3 = np.sqrt(C2_t3_x**2 + C2_t3_y**2 + C2_t3_z**2)


    recon_ds = xr.Dataset(data_vars={},
                   coords={
                              "x": ds.x.values,
                              "y": ds.y.values,
                              "z": ds.z.values,
                      },
                    attrs={"time":ds.time}
               )
    recon_ds["jx"] = jx
    recon_ds["jy"] = jy
    recon_ds["jz"] = jz
    recon_ds["mag_J"] = mag_J
    recon_ds["C2_t1_x"] = C2_t1_x
    recon_ds["C2_t1_y"] = C2_t1_y
    recon_ds["C2_t1_z"] = C2_t1_z
    recon_ds["mag_C2_t1"] = mag_C2_t1
    recon_ds["C2_t2_x"] = C2_t2_x
    recon_ds["C2_t2_y"] = C2_t2_y
    recon_ds["C2_t2_z"] = C2_t2_z
    recon_ds["mag_C2_t2"] = mag_C2_t2
    recon_ds["C2_t3_x"] = C2_t3_x
    recon_ds["C2_t3_y"] = C2_t3_y
    recon_ds["C2_t3_z"] = C2_t3_z
    recon_ds["mag_C2_t3"] = mag_C2_t3
    recon_ds.to_netcdf(f"data/{ds.time}_recon_metrics.nc")
    return


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

    time = str(pd.to_datetime(str(line[0]) + "_" + str(line[1]) + "_" + str(line[2]) + "_" + str(line[3]), format="%Y_%j_%H_%M"))

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
    ds.to_netcdf(f"data/{time}.nc")

    # If you wish to calculate the reconnection metrics, uncomment
    # the following function call

#    recon_metrics(ds)

    return

compute(omni_data[10])

# Set up process pool and operate the compute function on each entry
# in the omni_data array, or a subset of the array.
#with Pool(Nproc) as comp_pool:
#    comp_pool.map(compute,omni_data[0:8])
