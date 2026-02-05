from ctypes import CDLL, c_double, c_int
import ctypes as ct
import os
import numpy as np

libpath = os.path.dirname(os.path.realpath(__file__))
libfile = os.path.join(libpath, '../test.so')
ta16 = CDLL(libfile).ta16

ta16.argtypes = [ct.POINTER(c_int), ct.POINTER(c_double), ct.POINTER(c_double), ct.POINTER(c_double), ct.POINTER(c_double), ct.POINTER(c_double), ct.POINTER(c_double), ct.POINTER(c_double), c_int, c_int, c_int ]

date = np.array([2024, 133, 21, 0],dtype=c_int)

data = np.array([-829.0,57.8,91.8,0.2612,5.09,-87.7,1.9659,6.69],dtype=c_double)

datap = data.ctypes.data_as(ct.POINTER(ct.c_double))

datep = date.ctypes.data_as(ct.POINTER(ct.c_int))

x = np.linspace(-5,5,100)

y = np.linspace(-5,5,100)

z = np.linspace(-5,5,100)

bx = np.zeros( (len(x),len(y),len(z)), order= "F")

by = np.zeros( (len(x),len(y),len(z)), order= "F")

bz = np.zeros( (len(x),len(y),len(z)), order= "F")

xp = x.ctypes.data_as(ct.POINTER(ct.c_double))

yp = y.ctypes.data_as(ct.POINTER(ct.c_double))

zp = z.ctypes.data_as(ct.POINTER(ct.c_double))

bxp = bx.ctypes.data_as(ct.POINTER(ct.c_double))

byp = by.ctypes.data_as(ct.POINTER(ct.c_double))

bzp = bz.ctypes.data_as(ct.POINTER(ct.c_double))

print(bzp)

#ta16(datep,datap,xp,yp,zp,bxp,byp,bzp,ct.c_int(100),ct.c_int(100),ct.c_int(100))
