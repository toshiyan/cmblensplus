import libcurvedsky
import numpy

def qtt(est,lmax,rlmin,rlmax,TT,OCT,lfac=''):
  """
  Normalization of reconstructed fields from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TT [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_quad.qtt(est,lmax,rlmin,rlmax,TT,OCT,lfac):
  """
  return libcurvedsky.norm_quad.qtt(est,lmax,rlmin,rlmax,TT,OCT,lfac)

def qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,lfac=''):
  """
  Normalization of reconstructed fields from the TE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_quad.qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,lfac):
  """
  return libcurvedsky.norm_quad.qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,lfac)

def qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,lfac=''):
  """
  Normalization of reconstructed fields from the TB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :TE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_quad.qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,lfac):
  """
  return libcurvedsky.norm_quad.qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,lfac)

def qee(est,lmax,rlmin,rlmax,EE,OCE,lfac=''):
  """
  Normalization of reconstructed amplitude modulation from the EE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_quad.qee(est,lmax,rlmin,rlmax,EE,OCE,lfac):
  """
  return libcurvedsky.norm_quad.qee(est,lmax,rlmin,rlmax,EE,OCE,lfac)

def qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,lfac='',BB=0):
  """
  Normalization of reconstructed fields from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optionals): 
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_quad.qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,BB,lfac):
  """
  if BB==0: BB=0.*EE
  return libcurvedsky.norm_quad.qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,BB,lfac)

def qbb(est,lmax,rlmin,rlmax,BB,OCB,lfac=''):
  """
  Normalization of reconstructed amplitude modulation from the BB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :BB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Al [*2,l*] (*double*): Normalization, with bounds (0:lmax)

  Usage:
    :Al = curvedsky.norm_quad.qbb(est,lmax,rlmin,rlmax,BB,OCB,lfac):
  """
  return libcurvedsky.norm_quad.qbb(est,lmax,rlmin,rlmax,BB,OCB,lfac)

def qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,lfac=''):
  """
  Correlation between unnormalized TT and TE quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCTT [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :fCTE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = curvedsky.norm_quad.qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,lfac):
  """
  return libcurvedsky.norm_quad.qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,lfac)

def qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,lfac=''):
  """
  Correlation between unnormalized TT and EE quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCTT [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :fCEE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = curvedsky.norm_quad.qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,lfac):
  """
  return libcurvedsky.norm_quad.qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,lfac)

def qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,lfac=''):
  """
  Correlation between unnormalized TE and EE quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCEE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :fCTE [*l*] (*double*): Theory TE spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = curvedsky.norm_quad.qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,lfac):
  """
  return libcurvedsky.norm_quad.qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,lfac)

def qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,lfac=''):
  """
  Correlation between unnormalized TB and EB quadratic estimators

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCEE [*l*] (*double*): Theory EE spectrum, with bounds (0:rlmax)
    :fCBB [*l*] (*double*): Theory BB spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:rlmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)
    :OCTE [*l*] (*double*): Observed TE spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Ig [*l*] (*double*): Correlation between lensing potential estimators, with bounds (0:lmax)
    :Ic [*l*] (*double*): Correlation between curl mode estimators, with bounds (0:lmax)

  Usage:
    :Ig,Ic = curvedsky.norm_quad.qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,lfac):
  """
  return libcurvedsky.norm_quad.qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,lfac)

def qmv(lmax,QDO,Al,Il):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :lmax (*int*): Maximum multipole of the output power spectra
    :QDO[*6*] (*bool*): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. Currently, BB is always False
    :Al [*5,l*] (*double*): Normalizations of each estimator (TT, TE, EE, TB, EB).
    :Il [*4,l*] (*double*): Correlation between different estimators (TTxTE, TTxEE, TExEE, TBxEB).

  Returns:
    :MV [*l*] (*double*): Normalization of the MV estimator, with bounds (0:lmax)
    :Nl [*6,l*] (*double*): Weights for each estimator (TT, TE, EE, TB, EB, BB=0), with bounds (0:lmax)

  Usage:
    :MV,Nl = curvedsky.norm_quad.qmv(lmax,QDO,Al,Il):
  """
  return libcurvedsky.norm_quad.qmv(lmax,QDO,Al,Il)

def qall(est,QDO,lmax,rlmin,rlmax,fC,OC,lfac=''):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :QDO[*6*] (*bool*): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB.
    :lmax (*int*): Maximum multipole of the output power spectra
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC/OC [*l*] (*double*): Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Ag [*6,l*] (*double*): Normalization of the TT, TE, EE, TB, EB, and MV estimators for lensing potential, with bounds (6,0:lmax)
    :Ac [*6,l*] (*double*): Same as Ag but for curl mode
    :Nlg [*6,l*] (*double*): Weights for TT, TE, EE, TB, EB, and BB (=0) estimators for lensing potential, with bounds (6,0:lmax)
    :Nlc [*6,l*] (*double*): Same as Nlg but for curl mode

  Usage:
    :Ag,Ac,Nlg,Nlc = curvedsky.norm_quad.qall(est,QDO,lmax,rlmin,rlmax,fC,OC,lfac):
  """
  return libcurvedsky.norm_quad.qall(est,QDO,lmax,rlmin,rlmax,fC,OC,lfac)

def qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter=1,conv=0.00001):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization
    :elmax (*int*): Maximum multipole of input EE spectra, CE and OCE
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :dlmin/dlmax (*int*): Minimum/Maximum multipole of E mode and lensing potential for delensing
    :CE [*l*] (*double*): Theory EE angular power spectrum, with bounds (0:elmax)
    :OCE [*l*] (*double*): Observed EE spectrum, with bounds (0:elmax)
    :OCB [*l*] (*double*): Observed BB spectrum, with bounds (0:rlmax)
    :Cpp [*l*] (*double*): Theory lensing potential spectrum, with bounds (0:dlmax)

  Args(optional):
    :iter (*int*): number of iteration, default to 1 (no iteration)
    :conv (*double*): a parameter for convergence the iteration, default to 0.00001

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = curvedsky.norm_quad.qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter,conv):
  """
  return libcurvedsky.norm_quad.qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,iter,conv)

def xtt(est,lmax,rlmin,rlmax,fC,OCT,lfac=''):
  """
  Unnormalized response between lensing potential and amplitude modulation from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output normalization spectrum
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory TT spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed TT spectrum, with bounds (0:rlmax)

  Args(optional):
    :lfac (*str*): Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)

  Returns:
    :Ag [*l*] (*double*): Cross normalization, with bounds (0:lmax)

  Usage:
    :Ag = curvedsky.norm_quad.xtt(est,lmax,rlmin,rlmax,fC,OCT,lfac):
  """
  return libcurvedsky.norm_quad.xtt(est,lmax,rlmin,rlmax,fC,OCT,lfac)

def xeb(est,lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB=None):
  """
  Response of reconstructed field to other in the EB quadratic estimator

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
    :Aa = curvedsky.norm_quad.xeb(est,lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB):
  """
  if BB is None: BB= EE*0
  return libcurvedsky.norm_quad.xeb(est,lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB)

