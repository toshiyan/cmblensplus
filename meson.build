# Project definition.
# First argument is the project name. Second argument declares the languages used by the build.
# Here we need C because f2py generates C wrapper code, and Fortran because basic/flatsky are Fortran extensions.
project(
  'cmblensplus',
  ['c', 'fortran'],
  version: '0.1.0',
  default_options: [
    # Build optimized binaries by default.
    # For debugging, this can be changed to 'debug'.
    'buildtype=release',

    # Keep compiler warnings moderate.
    # Increase this if you want stricter checks.
    'warning_level=1',
  ],
)

# Load Meson's Python module.
# This provides helper functions for Python packages and extension modules.
py_mod = import('python')

# Find the Python interpreter used for this build.
# pure: false means this package can contain compiled extension modules.
py = py_mod.find_installation(pure: false)

# Load Meson's filesystem helper module.
# This is useful for path handling and file checks.
fs = import('fs')

# Absolute path to the source tree root.
# Use this for command arguments or external library paths,
# but do not pass source-tree absolute paths to include_directories().
project_root = meson.project_source_root()

# Absolute path to the build tree root.
# Meson writes generated files here, not into the source tree.
project_build_root = meson.project_build_root()

# Path to the script that converts .src files into .f90 files
# with f2py directives inserted.
#
# Use files() here instead of project_root / ...
# so Meson treats this as a source file in the project.
f2pysig = files('scripts/f2pysig.py')[0]

# Ask Python where NumPy's C headers are installed.
# This is needed to compile f2py-generated C wrapper code.
# This path is outside the project source tree, so it is OK as an absolute path.
numpy_inc = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()

# Ask Python where NumPy's f2py support sources are installed.
# This directory contains fortranobject.c and fortranobject.h,
# which are required when building f2py extension modules manually.
# This path is outside the project source tree, so it is OK as an absolute path.
numpy_f2py_src = run_command(
  py,
  [
    '-c',
    'import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))',
  ],
  check: true,
).stdout().strip()

# Create a Meson include directory object for NumPy and f2py headers.
# Subdirectories can reuse this variable when building libbasic/libflatsky.
#
# These are external include directories, not project source-tree directories.
numpy_dep = declare_dependency(
  compile_args: [
    '-I' + numpy_inc,
    '-I' + numpy_f2py_src,
  ],
)

# Get the Fortran compiler object.
# Subdirectories can use this to find Fortran libraries.
fc = meson.get_compiler('fortran')

# Path to the legacy/internal Fortran source root.
# This absolute path is fine for passing to find_library() dirs,
# but not for include_directories().
fortran_internal_root = project_root / 'fortran_internal'

# Directory where internal Fortran static/shared libraries are expected.
# For example:
#   fortran_internal/lib/libutils.a
#   fortran_internal/lib/libdft.a
#
# This is used by fc.find_library(..., dirs: [...]).
fortran_internal_lib_dir = project_root / 'fortran_internal' / 'lib'

# Directory containing the old Makefile for libutils.
# This Makefile builds Fortran utility routines used by libbasic.
fortran_internal_src_utils_dir = project_root / 'fortran_internal' / 'src_utils'

# Expected location of libutils.a after:
#   make -C fortran_internal/src_utils install
#
# This is the legacy install location used by the old build system.
fortran_internal_utils_lib = project_root / 'fortran_internal' / 'lib' / 'libutils.a'

# Fortran compiler selected by Meson.
# In your py_dev environment this is likely:
#   x86_64-conda-linux-gnu-gfortran
fc_cmd = fc.cmd_array()[0]

# Build libutils.a and install .mod files using the existing Makefile.
#
# This runs:
#   make -C fortran_internal/src_utils FC=<meson fortran compiler>
#   make -C fortran_internal/src_utils install FC=<meson fortran compiler>
#
# Expected side effects:
#   fortran_internal/mod/*.mod
#   fortran_internal/lib/libutils.a
#
# Meson formally tracks only the copied build-tree libutils.a.
utils_lib_target = custom_target(
  'fortran_internal_libutils',
  output: 'libutils.a',
  command: [
    'sh', '-c',
    'make -C "$1" FC="$4" && make -C "$1" install FC="$4" && cp "$2" "$3"',
    'sh',
    fortran_internal_src_utils_dir,
    fortran_internal_utils_lib,
    '@OUTPUT@',
    fc_cmd,
  ],
  install: false,
)

# Build-tree copy of libutils.a used for linking libbasic.
utils_lib_path = meson.current_build_dir() / 'libutils.a'


# Create an include directory object for internal Fortran module files.
#
# IMPORTANT:
# Do not use:
#   include_directories(fortran_internal_root / 'mod')
#
# because Meson rejects absolute paths into the source tree.
# Use a path relative to this meson.build file instead.
inc_fortran_internal = include_directories(
  'fortran_internal/mod',
)

# Install the top-level Python package file:
#   cmblensplus/__init__.py
#
# This ensures:
#   import cmblensplus
# works after installation.
py.install_sources(
  'cmblensplus/__init__.py',
  subdir: 'cmblensplus',
)

# Install the pure-Python curvedsky package.
# The file list is handled inside cmblensplus/curvedsky/meson.build.
subdir('cmblensplus/curvedsky')

# Build and install the f2py extension module for basic.
# The .src -> .f90 generation and libbasic build rules live there.
subdir('fortran_wrapped/src_basic')

# Build and install the f2py extension module for flatsky.
# The .src -> .f90 generation and libflatsky build rules live there.
subdir('fortran_wrapped/src_flatsky')