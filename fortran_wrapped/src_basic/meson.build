# Build the f2py extension module

# Name of the Python extension module.
basic_module_name = 'libbasic'

# Base names of .src files in this directory.
basic_src_names = [ 'wigner_funcs', 'bispec', 'flat', 'delens', 'galaxy', 'cosmofuncs', ]


# --------------------------------------------------------------------------- 
# Step 1: .src -> .f90 
# ---------------------------------------------------------------------------

# Convert each .src file into a generated .f90 file.
basic_generated_f90 = []

foreach name : basic_src_names
  basic_generated_f90 += custom_target(
    'basic_' + name + '_generated_f90',
    
    input: name + '.src', # input file
    output: name + '.f90', # output file

    # Run:
    #   python scripts/f2pysig.py -srcname <input> -o <output>
    command: [
      py,
      f2pysig,
      '-srcname', '@INPUT@',
      '-o', '@OUTPUT@',
    ],
    depends: [
      utils_lib_target,
    ],
  )
endforeach


# ---------------------------------------------------------------------------
# Step 2: generated .f90 -> f2py C/F90 wrapper sources
# ---------------------------------------------------------------------------

basic_f2py_wrappers = custom_target(
  'basic_f2py_wrappers',
  input: basic_generated_f90,
  output: [
    'basicmodule.c',
    'basic-f2pywrappers2.f90',
  ],
  command: [
    py,
    '-m', 'numpy.f2py',
    '--lower',
    '-m', basic_module_name,
    '@INPUT@',
    '--build-dir', '@OUTDIR@',
  ],
  depends: [
    utils_lib_target,
  ],
)


# ---------------------------------------------------------------------------
# Step 3: build Python extension module
# ---------------------------------------------------------------------------

py.extension_module(
  basic_module_name,
  # Sources for the extension:
  [
    basic_f2py_wrappers,
    basic_generated_f90,
    numpy_f2py_src / 'fortranobject.c',
  ],
  subdir: 'cmblensplus/basic/',

  include_directories: [
    inc_fortran_internal,
  ],

  dependencies: [
    numpy_dep,
  ],

  # Link directly to the libutils.a that Meson copied into the build dir.
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

