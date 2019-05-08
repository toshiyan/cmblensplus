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

