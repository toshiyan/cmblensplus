import libcurvedsky

def qtt(lmax,rlmin,rlmax,fC,OCT):
  """
  Normalization of reconstructed amplitude modulation from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Returns:
    :At [*l*] (*double*): tau normalization, with bounds (0:lmax)

  Usage:
    :At = curvedsky.norm_tau.qtt(lmax,rlmin,rlmax,fC,OCT):
  """
  return libcurvedsky.norm_tau.qtt(lmax,rlmin,rlmax,fC,OCT)

def qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB=None):
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
    :At = curvedsky.norm_tau.qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB):
  """
  if BB is None: BB= EE*0
  return libcurvedsky.norm_tau.qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB)

def oeb(lmax,rlmin,rlmax,EB,OCE,OCB):
  """
  Normalization of reconstructed amplitude from the EB quadratic estimator
  The kernels are the same as the rotation normalization but with a factor of 4. 

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

def stt(lmax,rlmin,rlmax,fC,OCT):
  """
  Unnormalized response between patchy tau and poisson sources/inhomogeneous noise with the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Returns:
    :At [*l*] (*double*): tau normalization, with bounds (0:lmax)

  Usage:
    :At = curvedsky.norm_tau.stt(lmax,rlmin,rlmax,fC,OCT):
  """
  return libcurvedsky.norm_tau.stt(lmax,rlmin,rlmax,fC,OCT)

