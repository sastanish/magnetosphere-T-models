import os
from numpy import ctypeslib, zeros
from ctypes import c_int, byref, util, CDLL

def ta16(date,data,x,y,z):
  nx = len(x)
  ny = len(y)
  nz = len(z)
  bx = zeros( (nx,ny,nz) )
  by = zeros( (nx,ny,nz) )
  bz = zeros( (nx,ny,nz) )

  ta16 = fortran_lib.ta16

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

if __name__=="__main__":
  lib_name = 'magnetomodels'

  if os.name == 'posix': # Linux or macOS
    if os.uname().sysname == 'Darwin': # macOS
      lib_path = util.find_library(lib_name) or os.path.join(os.path.expanduser('~/.local/lib'), f'lib{lib_name}.dylib')
    else: # Linux
      lib_path = util.find_library(lib_name) or os.path.join(os.path.expanduser('~/.local/lib'), f'lib{lib_name}.so')
  elif os.name == 'nt': # Windows
      lib_path = util.find_library(lib_name) or os.path.join(os.environ['LOCALAPPDATA'], 'Programs', 'Python', 'Python3X', 'Lib', 'site-packages', f'{lib_name}.dll')
  else:
    raise RuntimeError("Unsupported operating system")

  if not lib_path:
    raise FileNotFoundError(f"Could not find shared library '{lib_name}'. Make sure it's installed or in your system's library path.")

  try:
    fortran_lib = CDLL(lib_path)
  except OSError as e:
    print(f"Error loading library: {e}")
    print(f"Attempted path: {lib_path}")


