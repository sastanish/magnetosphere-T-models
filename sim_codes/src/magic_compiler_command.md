For the Intel compiler
```
FC="ifort" FFlags="qopenmp" python -m numpy.f2py -lgomp -c geopack.f reconnectionMetrics.f90 TA16_RBF.f TS04c.for TsyganenkoWrapper.f90 TsyganenkoWrapper.pyf -m TsyganenkoWrapper --backend meson

```
This is the command that works for the Gfort compiler.
```
FC="gfortran" FFlags="fopenmp" python -m numpy.f2py -lgomp -c geopack.f reconnectionMetrics.f90 TA16_RBF.f TS04c.for TsyganenkoWrapper.f90 TsyganenkoWrapper.pyf -m TsyganenkoWrapper --backend meson
```
