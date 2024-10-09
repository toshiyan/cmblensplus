import libcurvedsky
import numpy

def qte(est,lmax,rlmin,rlmax,TB,OCT,OCE,lfac=''):
  """
  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the TE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TB [*l*] (*double*): Theory TB spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_imag.qte(est,lmax,rlmin,rlmax,TB,OCT,OCE,lfac):
  """
  return libcurvedsky.norm_imag.qte(est,lmax,rlmin,rlmax,TB,OCT,OCE,lfac)

def qtb(est,lmax,rlmin,rlmax,TB,OCT,OCB,lfac=''):
  """
  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TB [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_imag.qtb(est,lmax,rlmin,rlmax,TB,OCT,OCB,lfac):
  """
  return libcurvedsky.norm_imag.qtb(est,lmax,rlmin,rlmax,TB,OCT,OCB,lfac)

def qee(est,lmax,rlmin,rlmax,fC,OCE,lfac=''):
  """
  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the E-mode quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_imag.qee(est,lmax,rlmin,rlmax,fC,OCE,lfac):
  """
  return libcurvedsky.norm_imag.qee(est,lmax,rlmin,rlmax,fC,OCE,lfac)

def qeb(est,lmax,rlmin,rlmax,fC,OCE,OCB,lfac=''):
  """
  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory EB spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_imag.qeb(est,lmax,rlmin,rlmax,fC,OCE,OCB,lfac):
  """
  return libcurvedsky.norm_imag.qeb(est,lmax,rlmin,rlmax,fC,OCE,OCB,lfac)

def qbb(est,lmax,rlmin,rlmax,fC,OCB,lfac=''):
  """
  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the B-mode quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_imag.qbb(est,lmax,rlmin,rlmax,fC,OCB,lfac):
  """
  return libcurvedsky.norm_imag.qbb(est,lmax,rlmin,rlmax,fC,OCB,lfac)

def qbb_asym(est,lmax,rlmin,rlmax,fC,OCB1,OCB2,lfac=''):
  """
  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the B-mode quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCB1/2 [*l*] (*double*): Observed BB spectrum for experiment 1 and 2, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_imag.qbb_asym(est,lmax,rlmin,rlmax,fC,OCB1,OCB2,lfac):
  """
  return libcurvedsky.norm_imag.qbb_asym(est,lmax,rlmin,rlmax,fC,OCB1,OCB2,lfac)

