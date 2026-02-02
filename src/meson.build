project('magnetosphere-T-models', 'fortran')

# Fortran main

# Python module
py = import('python').find_installation()

sources = files('core/geopack.f', 'core/inout.f90', 'core/TA16_RBF.f', 'core/TS05.f', 'core/control.f90', 'py_wrapper/python_wrapper.pyf', 'py_wrapper/python_wrapper.pyf')

add_mod = custom_target(
  'add_extension',
  input: sources,
  output: ['add' + py.extension_suffix()],
  command: [
    py.full_path(), '-m', 'numpy.f2py',
    '-c', 'core/geopack.f', 'core/inout.f90', 
    'core/TA16_RBF.f', 'core/TS05.f', 
    'core/control.f90', 'py_wrapper/python_wrapper.pyf', 
    'py_wrapper/python_wrapper.pyf', '-m', 'magnetosphere-T-models'
    ],
  build_by_default: true
)


# Install into build
install_subdir('.', install_dir: join_paths(meson.source_root(), 'builddir'),
               strip_directory: false,
               exclude_files: ['meson.build'])

# Also install the built extension (place it beside __init__.py)
install_data(add_mod, install_dir: join_paths(meson.source_root(), 'builddir'))
