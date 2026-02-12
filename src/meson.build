project('magnetomodels', 'fortran',
        version : '1.0',
        default_options : ['warning_level=1', 'python.install_env=auto']
       )

add_languages('fortran')
add_languages('c')

# Vars

netcdf_dep = dependency('netcdf', language: 'fortran')

omp_dep = dependency('openmp', required: false)

core_sources = files('core/geopack.f', 'core/inout.f90', 'core/TA16_RBF.f', 'core/TS05.f', 'core/control.f90', 'magnetomodels/ctype_wrapper.f90')
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

