import libcurvedsky

def qtb(lmax,rlmin,rlmax,fC,OCT,OCB):
  """
  Normalization of reconstructed pol. rot. angle from the TB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TE angular power spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed temperature angular power spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed B-mode angular power spectrum, with bounds (0:rlmax)

  Returns:
    :Aa [*l*] (*double*): Pol. rot. angle normalization, with bounds (0:lmax)

  Usage:
    :Aa = curvedsky.norm_rot.qtb(lmax,rlmin,rlmax,fC,OCT,OCB):
  """
  return libcurvedsky.norm_rot.qtb(lmax,rlmin,rlmax,fC,OCT,OCB)

def qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB=None):
  """
  Normalization of reconstructed pol. rot. angle from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optionals): 
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)

  Returns:
    :Aa [*l*] (*double*): Pol. rot. angle normalization, with bounds (0:lmax)

  Usage:
    :Aa = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB):
  """
  if BB is None: BB= EE*0
  return libcurvedsky.norm_rot.qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB)

def teb(lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB=None):
  """
  Response of reconstructed pol. rot. angle to amplitude in the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :EB [*l*] (*double*): Theory EB spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optionals): 
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)

  Returns:
    :Aa [*l*] (*double*): Pol. rot. angle normalization, with bounds (0:lmax)

  Usage:
    :Aa = curvedsky.norm_rot.teb(lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB):
  """
  if BB is None: BB= EE*0
  return libcurvedsky.norm_rot.teb(lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB)

