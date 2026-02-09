project('magnetosphere-T-models', 'fortran',
        version : '1.0',
        default_options : ['warning_level=1']
       )

add_languages('fortran')
add_languages('c')

# Vars

netcdf_dep = dependency('netcdf', language: 'fortran')

omp_dep = dependency('openmp', required: false)

core_sources = files('core/geopack.f', 'core/inout.f90', 'core/TA16_RBF.f', 'core/TS05.f', 'core/control.f90', 'python_wrapper/ctype_wrapper.f90')
core = library('magnetomodels', 
               core_sources, 
               dependencies : [netcdf_dep, omp_dep],
               install : true)

# Fortran main
executable('main_ta16.out',
           'main_program.f90',
           link_with: core,
           dependencies : [netcdf_dep, omp_dep]
           )

# Python module

#py_mod = import('python')
#py = py_mod.find_installation(pure: false)
#py_dep = py.dependency()

#incdir_numpy = run_command(py,
#  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
#  check : true
#).stdout().strip()


#py.extension_module(
#  'magnetoTmodels',
#  'py_wrapper/_magnetoTmodels.c',
#  dependencies : py_dep,
#  include_directories : incdir_numpy,
#  link_with: [core],#, lib_f90_wrap],
#  install : true
#)
