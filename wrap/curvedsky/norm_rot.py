import libcurvedsky

def qeb(lmax,rlmin,rlmax,fC,OCE,OCB):
  """
  Normalization of reconstructed pol. rot. angle from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory EE angular power spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed E-mode angular power spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed B-mode angular power spectrum, with bounds (0:rlmax)

  Returns:
    :Aa [*l*] (*double*): Pol. rot. angle normalization, with bounds (0:lmax)

  Usage:
    :Aa = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,fC,OCE,OCB):
  """
  return libcurvedsky.norm_rot.qeb(lmax,rlmin,rlmax,fC,OCE,OCB)

