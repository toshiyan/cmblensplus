import libcurvedsky

def qtt(lmax,rlmin,rlmax,OCT):
  """
  Normalization of reconstructed src field from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory temperature angular power spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed temperature angular power spectrum, with bounds (0:rlmax)

  Returns:
    :As [*l*] (*double*): src field normalization, with bounds (0:lmax)

  Usage:
    :As = curvedsky.norm_src.qtt(lmax,rlmin,rlmax,OCT):
  """
  return libcurvedsky.norm_src.qtt(lmax,rlmin,rlmax,OCT)

