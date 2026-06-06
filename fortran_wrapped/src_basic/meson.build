# Build the f2py extension module:
#   libbasic
#
# This file is read from the top-level meson.build via:
#   subdir('fortran_wrapped/src_basic')
#
# Variables assumed to be defined in the top-level meson.build:
#   py
#   fc
#   f2pysig
#   numpy_f2py_src
#   inc_numpy
#   inc_fortran_internal
#   fortran_internal_lib_dir


# Name of the Python extension module.
# This preserves the old import style:
#   import libbasic
basic_module_name = 'libbasic'


# Base names of .src files in this directory.
# Old setup.py used:
#   wigner_funcs, bispec, flat, delens, galaxy, cosmofuncs
#
# Each name corresponds to:
#   <name>.src  -> generated <name>.f90
basic_src_names = [
  'wigner_funcs',
  'bispec',
  'flat',
  'delens',
  'galaxy',
  'cosmofuncs',
]


# Convert each .src file into a generated .f90 file.
#
# The generated .f90 files are written into Meson's build directory,
# not into the source tree.
basic_generated_f90 = []

foreach name : basic_src_names
  basic_generated_f90 += custom_target(
    'basic_' + name + '_generated_f90',

    # Input file, for example:
    #   wigner_funcs.src
    input: name + '.src',

    # Output file in the build directory, for example:
    #   wigner_funcs.f90
    output: name + '.f90',

    # Run:
    #   python scripts/f2pysig.py -srcname <input> -o <output>
    command: [
      py,
      f2pysig,
      '-srcname', '@INPUT@',
      '-o', '@OUTPUT@',
    ],
  )
endforeach



# Build the Python extension module.
#
# This creates something like:
#   libbasic.cpython-313-x86_64-linux-gnu.so
#
# It can be imported as:
#   import libbasic
py.extension_module(
  basic_module_name,

  # Sources for the extension:
  #   1. generated Fortran files from .src
  #   2. NumPy f2py runtime C source
  [
    basic_generated_f90,
    numpy_f2py_src / 'fortranobject.c',
  ],

  include_directories: [
    inc_fortran_internal,
  ],

  dependencies: [
    numpy_dep,
  ],

  # Link directly to the libutils.a that Meson copied into the build dir.
  #
  # We do not use fc.find_library('utils') because libutils.a does not exist
  # at Meson configure time. It is created later by utils_lib_target.
  link_args: [
    utils_lib_path,
  ],

  # Ensure libutils.a is built before linking libbasic.
  link_depends: [
    utils_lib_target,
  ],

  # Optimization flags for Fortran compilation.
  fortran_args: [
    '-O3',
  ],

  # Install the extension into site-packages.
  install: true,
)

