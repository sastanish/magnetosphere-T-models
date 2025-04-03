# 2024 may storm
# Read OMNI data
# Calculate magnetic field with TS04
# Determine reconnection quantities by finite differences
# Output xarray_dataset file for visualization

import numpy as np
import pandas as pd
import xarray as xr
import TS04c
import geopack_dipoles
from alive_progress import alive_bar

# Part 1: read data
# The data has been pre-selected and placed in a file called TS04_2017_3to13_sep

omni_data = np.genfromtxt('may_2024_storm_with_TS05_vars.dat',dtype=None)

# Setup the desired GSW Coordinates and data-structure

(nx, ny, nz) = (60, 200, 200)
x0 = -10
x1 = 2
y0 = -10
y1 = 10
z0 = -10
z1 = 10

x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)
z = np.linspace(z0, z1, nz)

times = [ str(d[0]) + "_" + str(d[1]) + "_" + str(d[2]) + "_" + str(d[3])
         for d in omni_data[0:-1:12]]
dtime = pd.to_datetime(times, format="%Y_%j_%H_%M")

[bx, by, bz] = np.zeros( (3, len(x), len(y), len(z), len(times)) )

# Part 2: Calculate the field via TS04

with alive_bar(bx.size) as bar:
    for i,line in enumerate(omni_data[0:-1:12]):
        iopt = 0
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
        for ix, val_x in enumerate(x):
            for iy, val_y in enumerate(y):
                for iz, val_z in enumerate(z):
                    [extern_bx, extern_by, extern_bz] = TS04c.t04_s(iopt,parmod,ps,val_x,val_y,val_z)
                    [intern_bx, intern_by, intern_bz] = geopack_dipoles.igrf_gsw_08(iyear,iday,ihour,imin,isec,vgsex,vgsey,vgsez,val_x,val_y,val_z)
#                    [intern_bx, intern_by, intern_bz] = geopack_dipoles.dip_08(iyear,iday,ihour,imin,isec,vgsex,vgsey,vgsez,val_x,val_y,val_z)
                    bx[ix,iy,iz,i] = extern_bx + intern_bx
                    by[ix,iy,iz,i] = extern_by + intern_by
                    bz[ix,iy,iz,i] = extern_bz + intern_bz
                    bar()
# Still getting a vtk writer to work
#        write_vtk( (bx[...,i],by[...,i],bz[...,i]), dtime[i], file_name="vtk_files/SeptemberStorm")


dataset = xr.Dataset(
                data_vars={
                            "bx": (["x", "y", "z", "time"], bx),
                            "by": (["x", "y", "z", "time"], by),
                            "bz": (["x", "y", "z", "time"], bz),
                            },
                coords={
                                "x": x,
                                "y": y,
                                "z": z,
                        "omni_time": times,
                             "time": dtime,
                        }
                     )


# Write to several file formats

## netCDF
dataset.to_netcdf("MayStorm_2024_1h_cadence_igrf_08_dipole_magnetosheath.nc")

## VTK writer

def write_vtk(vecs,time,file_name="vtk_files/data"):
    [bx,by,bz] = vecs

    dim=bx.shape

    # Generate the grid
    xx,yy,zz=np.mgrid[0:dim[0],0:dim[1],0:dim[2]]
    pts = np.empty(dim + (3,), dtype=int)
    pts[..., 0] = xx
    pts[..., 1] = yy
    pts[..., 2] = zz

    vectors = np.empty(dim + (3,), dtype=float)
    vectors[..., 0] = bx
    vectors[..., 1] = by
    vectors[..., 2] = bz

    # We reorder the points and vectors so this is as per VTK's
    # requirement of x first, y next and z last.
    pts = pts.transpose(2, 1, 0, 3).copy()
    pts.shape = pts.size // 3, 3

    vectors = vectors.transpose(2, 1, 0, 3).copy()
    vectors.shape = vectors.size // 3, 3

    sg = vtkStructuredGrid(dimensions=xx.shape, points=pts)

    sg.point_data.vectors = vectors
    sg.point_data.vectors.name = 'B'

    write_data(sg, data_dir + str(time) + '.vtk')

