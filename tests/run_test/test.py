import ctypes as ct
import numpy as np
import os

os.system('export OMP_NUM_THREADS=1')

api = np.ctypeslib.load_library('libmagnetomodels.so','/usr/local/lib')

nx = 50
ny = 50
nz = 50

x = np.linspace(-5,5,nx)
y = np.linspace(-5,5,ny)
z = np.linspace(-5,5,nz)

bx = np.zeros( (nx,ny,nz) )
by = np.zeros( (nx,ny,nz) )
bz = np.zeros( (nx,ny,nz) )

omni_data = np.genfromtxt('./input_data.lst',dtype=None)

print(omni_data.shape)

date = np.array([omni_data[1][i] for i in range(4)])
data = np.array([omni_data[1][i] for i in [7,8,9,17,18,-1,-2,-4] ])

api.ta16( ct.byref(np.ctypeslib.as_ctypes(date)),
          ct.byref(np.ctypeslib.as_ctypes(data)),
          ct.byref(np.ctypeslib.as_ctypes(x)),
          ct.byref(np.ctypeslib.as_ctypes(y)),
          ct.byref(np.ctypeslib.as_ctypes(z)),
          ct.byref(np.ctypeslib.as_ctypes(bx)),
          ct.byref(np.ctypeslib.as_ctypes(by)),
          ct.byref(np.ctypeslib.as_ctypes(bz)),
          ct.byref(ct.c_int(nx)),
          ct.byref(ct.c_int(ny)),
          ct.byref(ct.c_int(nz))
          )
