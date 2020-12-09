import libcurvedsky

def norm_qtt(lmax,rlmin,rlmax,TT,OCT,est):
  """
  Usage:
    :Al = curvedsky.norm_general.norm_qtt(lmax,rlmin,rlmax,TT,OCT,est):
  """
  return libcurvedsky.norm_general.norm_qtt(lmax,rlmin,rlmax,TT,OCT,est)

def norm_qte(lmax,rlmin,rlmax,TE,OCT,OCE,est):
  """
  Normalization of reconstructed amplitude modulation from the TE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Returns:
    :Al [*l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_general.norm_qte(lmax,rlmin,rlmax,TE,OCT,OCE,est):
  """
  return libcurvedsky.norm_general.norm_qte(lmax,rlmin,rlmax,TE,OCT,OCE,est)

def norm_qtb(lmax,rlmin,rlmax,TB,OCT,OCB):
  """
  Normalization of reconstructed amplitude modulation from the TB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TB [*l*] (*double*): Theory TB spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Returns:
    :Al [*l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_general.norm_qtb(lmax,rlmin,rlmax,TB,OCT,OCB):
  """
  return libcurvedsky.norm_general.norm_qtb(lmax,rlmin,rlmax,TB,OCT,OCB)

def norm_qee(lmax,rlmin,rlmax,EE,OCE):
  """
  Normalization of reconstructed amplitude modulation from the EE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Returns:
    :Al [*l*] (*double*): Amplitude modulation normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_general.norm_qee(lmax,rlmin,rlmax,EE,OCE):
  """
  return libcurvedsky.norm_general.norm_qee(lmax,rlmin,rlmax,EE,OCE)

def norm_qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB=None):
  """
  Normalization of reconstructed amplitude modulation from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optionals): 
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)

  Returns:
    :At [*l*] (*double*): Amplitude modulation normalization, with bounds (0:lmax)

  Usage:
    :At = curvedsky.norm_general.norm_qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB):
  """
  if BB is None: BB= EE*0
  return libcurvedsky.norm_general.norm_qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB)

def norm_qbb(lmax,rlmin,rlmax,BB,OCB):
  """
  Normalization of reconstructed amplitude modulation from the BB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Returns:
    :Al [*l*] (*double*): Amplitude modulation normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_general.norm_qbb(lmax,rlmin,rlmax,BB,OCB):
  """
  return libcurvedsky.norm_general.norm_qbb(lmax,rlmin,rlmax,BB,OCB)

