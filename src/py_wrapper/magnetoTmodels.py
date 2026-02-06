from __future__ import print_function, absolute_import, division
import _magnetoTmodels
import f90wrap.runtime
import logging
import numpy
import warnings
import weakref

class Compute(f90wrap.runtime.FortranModule):
    """
    Module compute
    Defined at python_wrapper.f90 lines 1-73
    """
    @staticmethod
    def ta16(date, data, x, y, z, bx, by, bz, nx, ny, nz, interface_call=False):
        """
        ta16(date, data, x, y, z, bx, by, bz, nx, ny, nz)
        Defined at python_wrapper.f90 lines 4-29
        
        Parameters
        ----------
        date : int array
        data : float array
        x : float array
        y : float array
        z : float array
        bx : float array
        by : float array
        bz : float array
        nx : int32
        ny : int32
        nz : int32
        """
        _magnetoTmodels.f90wrap_compute__ta16(date=date, data=data, x=x, y=y, z=z, \
            bx=bx, by=by, bz=bz, nx=nx, ny=ny, nz=nz)
    
    _dt_array_initialisers = []
    

compute = Compute()

