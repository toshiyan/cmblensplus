# Install pure-Python curvedsky modules.
# These files are copied into:
#   <site-packages>/cmblensplus/curvedsky/
#
# curvedsky no longer builds a Fortran extension in the new structure.

py.install_sources(
  # Package initializer.
  # Enables `import cmblensplus.curvedsky`.
  '__init__.py',

  # Public API modules.
  # Keep these names stable for backward compatibility.
  'rec_lens.py',
  'norm_quad.py',
  'norm_imag.py',
  'rec_tau.py',
  'rec_rot.py',
  'rec_src.py',
  'rec_iamp.py',
  'rec_ilens.py',
  'delens.py',
  'bispec.py',
  'cninv.py',
  'utils.py',

  # Install destination inside the Python package.
  subdir: 'cmblensplus/curvedsky',
)
