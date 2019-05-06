import libcurvedsky

def qtt(lmax,rlmin,rlmax,fC,OCT,gtype=None):
  """
  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Theory temperature angular power spectrum, with bounds (0:rlmax)
    :OCT [*l*] (*double*): Observed temperature angular power spectrum, with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag [*l*] (*double*): CMB lensing potential normalization, with bounds (0:lmax)
    :Ac [*l*] (*double*): Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)

  Usage:
    :Ag,Ac = curvedsky.norm_lens.qtt(lmax,rlmin,rlmax,fC,OCT,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qtt(lmax,rlmin,rlmax,fC,OCT,gtype)

def qte(lmax,rlmin,rlmax,fC,OCT,OCE,gtype=None):
  """
  Usage:
    :Ag,Ac = curvedsky.norm_lens.qte(lmax,rlmin,rlmax,fC,OCT,OCE,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qte(lmax,rlmin,rlmax,fC,OCT,OCE,gtype)

def qtb(lmax,rlmin,rlmax,fC,OCT,OCB,gtype=None):
  """
  Usage:
    :Ag,Ac = curvedsky.norm_lens.qtb(lmax,rlmin,rlmax,fC,OCT,OCB,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qtb(lmax,rlmin,rlmax,fC,OCT,OCB,gtype)

def qee(lmax,rlmin,rlmax,fC,OCE,gtype=None):
  """
  Usage:
    :Ag,Ac = curvedsky.norm_lens.qee(lmax,rlmin,rlmax,fC,OCE,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qee(lmax,rlmin,rlmax,fC,OCE,gtype)

def qeb(lmax,rlmin,rlmax,fC,OCE,OCB,gtype=None):
  """
  Usage:
    :Ag,Ac = curvedsky.norm_lens.qeb(lmax,rlmin,rlmax,fC,OCE,OCB,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qeb(lmax,rlmin,rlmax,fC,OCE,OCB,gtype)

def qbb(lmax,rlmin,rlmax,fC,OCB,gtype=None):
  """
  Usage:
    :Ag,Ac = curvedsky.norm_lens.qbb(lmax,rlmin,rlmax,fC,OCB,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qbb(lmax,rlmin,rlmax,fC,OCB,gtype)

def qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype=None):
  """
  Usage:
    :Ig,Ic = curvedsky.norm_lens.qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,gtype)

def qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype=None):
  """
  Usage:
    :Ig,Ic = curvedsky.norm_lens.qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,gtype)

def qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype=None):
  """
  Usage:
    :Ig,Ic = curvedsky.norm_lens.qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,gtype)

def qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype=None):
  """
  Usage:
    :Ig,Ic = curvedsky.norm_lens.qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,gtype)

def qmv(lmax,QDO,Al,Il):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :lmax (*int*): Maximum multipole of the output power spectra
    :QDO[*6*] (*bool*): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB.
    :Al [*5,l*] (*double*): Normalizations of each estimator (TT, TE, EE, TB, EB).
    :Il [*4,l*] (*double*): Correlation between different estimators (TTxTE, TTxEE, TExEE, TBxEB).

  Returns:
    :MV [*l*] (*double*): Normalization of the MV estimator, with bounds (0:lmax)
    :Nl [*6,l*] (*double*): Weights for each estimator (TT, TE, EE, TB, EB, BB=0), with bounds (0:lmax)

  Usage:
    :MV,Nl = curvedsky.norm_lens.qmv(lmax,QDO,Al,Il):
  """
  return libcurvedsky.norm_lens.qmv(lmax,QDO,Al,Il)

def qall(QDO,lmax,rlmin,rlmax,fC,OC,gtype=None):
  """
  Compute MV estimator normalization. Currently BB is ignored. 

  Args:
    :QDO[*6*] (*bool*): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB.
    :lmax (*int*): Maximum multipole of the output power spectra
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC/OC [*l*] (*double*): Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :Ag/Ac [*6,l*] (*double*): Normalization of the TT, TE, EE, TB, EB, and MV estimators for gradient/curl modes, with bounds (6,0:lmax)
    :Nlg/Nlc [*6,l*] (*double*): Weights for TT, TE, EE, TB, EB, and BB (=0) estimators for gradient/curl modes, with bounds (6,0:lmax)

  Usage:
    :Ag,Ac,Nlg,Nlc = curvedsky.norm_lens.qall(QDO,lmax,rlmin,rlmax,fC,OC,gtype):
  """
  if gtype is None: gtype= ''
  return libcurvedsky.norm_lens.qall(QDO,lmax,rlmin,rlmax,fC,OC,gtype)

def qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter=None,conv=None):
  """
  Usage:
    :Alg,Alc = curvedsky.norm_lens.qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter,conv):
  """
  if iter is None: iter= 1
  if conv is None: conv= 0.001
  return libcurvedsky.norm_lens.qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter,conv)

