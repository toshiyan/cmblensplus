
project(
  'cmblensplus',
  ['c', 'fortran'],
  version: '0.1.0',
  default_options: [
    'buildtype=release',
    'warning_level=1',
  ],
)

# ---------------------------------------------------------------------------
# Python / NumPy / compiler setup
# ---------------------------------------------------------------------------

py_mod = import('python')
py = py_mod.find_installation(pure: false)

project_root = meson.project_source_root()
fc = meson.get_compiler('fortran')
fc_cmd = fc.cmd_array()[0]

# Script used to convert .src files into .f90 files with f2py directives.
f2pysig = files('scripts/f2pysig.py')[0]

# NumPy headers and f2py support sources.
numpy_inc = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()

numpy_f2py_src = run_command(
  py,
  [
    '-c',
    'import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))',
  ],
  check: true,
).stdout().strip()

numpy_dep = declare_dependency(
  compile_args: [
    '-I' + numpy_inc,
    '-I' + numpy_f2py_src,
  ],
)

# Internal Fortran module directory.
#
# The legacy Makefiles install .mod files into:
#   fortran_internal/mod
inc_fortran_internal = include_directories(
  'fortran_internal/mod',
)

# ---------------------------------------------------------------------------
# External dependencies
# ---------------------------------------------------------------------------

# FFTW is required by fortran_internal/src_dft and flatsky wrappers.
#
# In conda/HPC/system environments, this is usually found via pkg-config.
# If pkg-config metadata is unavailable, fall back to direct library lookup.
fftw_dep = dependency('fftw3', required: false)

if not fftw_dep.found()
  fftw_dep = fc.find_library('fftw3', required: true)
endif

# FFTW compile flags for the legacy src_dft Makefile.
#
# dependency('fftw3') passes include flags to Meson targets, but not to external
# Makefile invocations. Therefore we explicitly query pkg-config when available.
fftw_pkg_cflags = run_command(
  'pkg-config',
  '--cflags',
  'fftw3',
  check: false,
).stdout().strip()

# ---------------------------------------------------------------------------
# Legacy internal Fortran libraries
# ---------------------------------------------------------------------------
#
# The old build system builds internal Fortran libraries with Makefiles:
#
#   fortran_internal/src_utils -> libutils.a + constants.mod, general.mod, ...
#   fortran_internal/src_dft   -> libdft.a   + fftw.mod, dft-related modules
#
# The dft library depends on modules produced by utils, e.g. constants.mod.
# Therefore the build order must be:
#
#   libutils.a first
#   libdft.a second
#   f2py extension modules last
#
# Meson formally tracks the copied build-tree .a files as custom_target outputs.
# The .mod files are side effects installed into fortran_internal/mod by the
# legacy Makefiles.

fortran_internal_lib_dir = project_root / 'fortran_internal' / 'lib'

# ---- libutils.a ------------------------------------------------------------

fortran_internal_src_utils_dir = project_root / 'fortran_internal' / 'src_utils'
fortran_internal_utils_lib = fortran_internal_lib_dir / 'libutils.a'

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

utils_lib_path = meson.current_build_dir() / 'libutils.a'

# ---- libdft.a --------------------------------------------------------------
#
# src_dft/fftw.f90 uses modules from src_utils.
# Therefore dft_lib_target must depend on utils_lib_target.

fortran_internal_src_dft_dir = project_root / 'fortran_internal' / 'src_dft'
fortran_internal_dft_lib = fortran_internal_lib_dir / 'libdft.a'

dft_fflags = '-O3 -cpp -fPIC -I../mod ' + fftw_pkg_cflags

dft_lib_target = custom_target(
  'fortran_internal_libdft',
  output: 'libdft.a',
  command: [
    'sh', '-c',
    'srcdir="$1"; libsrc="$2"; libout="$3"; fc="$4"; fflags="$5"; make -C "$srcdir" FC="$fc" FFLAGS="$fflags" && make -C "$srcdir" install FC="$fc" FFLAGS="$fflags" && cp "$libsrc" "$libout"',
    'sh',
    fortran_internal_src_dft_dir,
    fortran_internal_dft_lib,
    '@OUTPUT@',
    fc_cmd,
    dft_fflags,
  ],
  depends: [
    utils_lib_target,
  ],
  install: false,
)

dft_lib_path = meson.current_build_dir() / 'libdft.a'

# ---------------------------------------------------------------------------
# Python package installation
# ---------------------------------------------------------------------------

py.install_sources(
  'cmblensplus/__init__.py',
  subdir: 'cmblensplus',
)

subdir('cmblensplus/basic')
subdir('cmblensplus/curvedsky')
subdir('cmblensplus/flatsky')


# ---------------------------------------------------------------------------
# f2py extension modules
# ---------------------------------------------------------------------------
#
# These subdirectories use the variables defined above:
#
#   py
#   f2pysig
#   numpy_dep
#   numpy_f2py_src
#   inc_fortran_internal
#   utils_lib_target
#   utils_lib_path
#   dft_lib_target
#   dft_lib_path
#   fftw_dep
#
# In particular:
#
#   src_basic   should depend on utils_lib_target.
#   src_flatsky should depend on both utils_lib_target and dft_lib_target.

subdir('fortran_wrapped/src_basic')
subdir('fortran_wrapped/src_flatsky')


