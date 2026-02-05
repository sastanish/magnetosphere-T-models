project('magnetosphere-T-models', 'fortran',
        version : '1.0',
        default_options : ['warning_level=1']
       )

# Vars

netcdf_dep = dependency('netcdf', language: 'fortran')

omp_dep = dependency('openmp', required: false)

fortran_sources = files('core/geopack.f', 'core/inout.f90', 'core/TA16_RBF.f', 'core/TS05.f', 'core/control.f90', 'main_program.f90')

# Fortran main
add_languages('fortran')
executable('main_ta16.out',
           fortran_sources,
           dependencies : [netcdf_dep, omp_dep])

# Python module
#py_mod = import('python')
#py = py_mod.find_installation(pure: false)
#py_dep = py.dependency()

#incdir_numpy = run_command(py,
#                 ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
#                 check : true).stdout().strip()

#incdir_f2py = run_command(py,
#                ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
#                check : true).stdout().strip()

#incdir_netcdf = run_command(

#inc_cp = include_directories(incdir_numpy, incdir_f2py)

#magneto_source = custom_target(
#  'magnetosphere-T-module.c',
#  input : ['core/geopack.f', 'core/inout.f90', 'core/TA16_RBF.f', 'core/TS05.f', 'core/control.f90', 'py_wrapper/python_wrapper.pyf', 'py_wrapper/python_wrapper.f90'],
#  output : ['magnetosphere-T-module.c', 'magnetosphere-T-f2pywrappers.f'],
#  command : [py, '-m', 'numpy.f2py', '-c', '@INPUT@', '-m', 'magnetosphere-T-models', '--lower']
#  )

#py.extension_module(
#  'magnetosphere-T-models',
#  magneto_source,
#  incdir_f2py / 'fortranobject.c',
#  dependencies : py_dep,
#  install : true
#)
