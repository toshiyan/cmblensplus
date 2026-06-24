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


# Install private pure-Python backend modules replacing libcurvedsky.
py.install_sources(
  '_core/__init__.py',
  '_core/lib_bispec.py',
  '_core/lib_cninv.py',
  '_core/lib_delens.py',
  '_core/lib_rec_iamp.py',
  '_core/lib_rec_ilens.py',
  '_core/lib_rec_lens.py',
  '_core/lib_rec_rot.py',
  '_core/lib_rec_src.py',
  '_core/lib_rec_tau.py',
  '_core/lib_utils.py',
  '_core/cgd.py',
  '_core/sht_ducc0.py',
  subdir: 'cmblensplus/curvedsky/_core',
)
