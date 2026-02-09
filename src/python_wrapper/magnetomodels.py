import os
from numpy import ctypeslib, zeros
from ctypes import c_int, byref

def ta16(date,data,x,y,z):
  nx = len(x)
  ny = len(y)
  nz = len(z)
  bx = zeros( (nx,ny,nz) )
  by = zeros( (nx,ny,nz) )
  bz = zeros( (nx,ny,nz) )

  ta16 = ctypeslib.load_library('libmagnetomodels.so','/usr/local/lib')

  api.ta16( byref(ctypeslib.as_ctypes(date)),
            byref(ctypeslib.as_ctypes(data)),
            byref(ctypeslib.as_ctypes(x)),
            byref(ctypeslib.as_ctypes(y)),
            byref(ctypeslib.as_ctypes(z)),
            byref(ctypeslib.as_ctypes(bx)),
            byref(ctypeslib.as_ctypes(by)),
            byref(ctypeslib.as_ctypes(bz)),
            byref(c_int(nx)),
            byref(c_int(ny)),
            byref(c_int(nz))
            )
  return (bx, by, bz)
