project('magnetosphere-T-models', 'fortran',
        version : '0.1',
        default_options : ['warning_level=1']
       )

add_languages('fortran')

# Vars
netcdf_dep = dependency('netcdf', language: 'fortran')
omp_dep = dependency('openmp', required: false)
lib_install_dir = join_paths(meson.global_source_root(), 'lib')
bin_install_dir = join_paths(meson.global_source_root(), 'bin')

# Main Fortran Library
core_files = files('src/core/geopack.f', 'src/core/inout.f90', 'src/core/TA16_RBF.f', 'src/core/TS05.f', 'src/core/control.f90')
core = static_library('magnetosphere-T-models', 
               core_files, 
               dependencies : [netcdf_dep, omp_dep],
               install : true,
               install_dir : lib_install_dir)

# Fortran Main Files
executable('run_ta16.out',
           'src/main/run_ta16.f90',
           link_with: core,
           dependencies : [netcdf_dep, omp_dep],
           install : true,
           install_dir : bin_install_dir
           )
executable('run_ts05.out',
           'src/main/run_ts05.f90',
           link_with: core,
           dependencies : [netcdf_dep, omp_dep],
           install : true,
           install_dir : bin_install_dir
           )

# Python Wrapper

py_wrapper_files = files('src/core/geopack.f', 'src/core/TA16_RBF.f', 'src/core/TS05.f', 'src/core/control.f90', 'src/python/magnetoTpython.f90')
py_wrapper = static_library('magnetoTpython',
                            py_wrapper_files, 
                            dependencies : [omp_dep],
                            install : true,
                            install_dir : lib_install_dir)

# Python Library
py_mod = import('python')
py = py_mod.find_installation(pure : false)

py.install_sources(['src/python/magnetoTmodels.py'],
                   install_dir : join_paths(bin_install_dir,'python')
                   )

