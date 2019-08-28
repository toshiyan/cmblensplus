import libcurvedsky

def qtt(lmax,rlmin,rlmax,fC,OCT):
  """
  Normalization of reconstructed tau from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory temperature angular power spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed temperature angular power spectrum, with bounds (0:rlmax)

  Returns:
    :At [*l*] (*double*): tau normalization, with bounds (0:lmax)

  Usage:
    :At = curvedsky.norm_tau.qtt(lmax,rlmin,rlmax,fC,OCT):
  """
  return libcurvedsky.norm_tau.qtt(lmax,rlmin,rlmax,fC,OCT)

def oeb(lmax,rlmin,rlmax,EB,OCE,OCB):
  """
  Normalization of reconstructed amplitude from the EB quadratic estimator
  The kernels are the same as the rotation normalization but with a factor 4. 

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EB [*l*] (*double*): Theory EB spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Returns:
    :At [*l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :At = curvedsky.norm_tau.oeb(lmax,rlmin,rlmax,EB,OCE,OCB):
  """
  return libcurvedsky.norm_tau.oeb(lmax,rlmin,rlmax,EB,OCE,OCB)

