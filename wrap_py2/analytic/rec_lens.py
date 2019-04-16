import libanalytic

def qtt(lmax,rlmin,rlmax,fC,OCT):
  """
  Usage:
    - e.g., Ag,Ac = analytic.rec_lens.qtt(lmax,rlmin,rlmax,fC,OCT)
  """
  return analytic.rec_lens.qtt(lmax,rlmin,rlmax,fC,OCT)

def qte(lmax,rlmin,rlmax,fC,OCT,OCE):
  """
  Usage:
    - e.g., Ag,Ac = analytic.rec_lens.qte(lmax,rlmin,rlmax,fC,OCT,OCE)
  """
  return analytic.rec_lens.qte(lmax,rlmin,rlmax,fC,OCT,OCE)

def qtb(lmax,rlmin,rlmax,fC,OCT,OCB):
  """
  Usage:
    - e.g., Ag,Ac = analytic.rec_lens.qtb(lmax,rlmin,rlmax,fC,OCT,OCB)
  """
  return analytic.rec_lens.qtb(lmax,rlmin,rlmax,fC,OCT,OCB)

def qee(lmax,rlmin,rlmax,fC,OCE):
  """
  Usage:
    - e.g., Ag,Ac = analytic.rec_lens.qee(lmax,rlmin,rlmax,fC,OCE)
  """
  return analytic.rec_lens.qee(lmax,rlmin,rlmax,fC,OCE)

def qeb(lmax,rlmin,rlmax,fC,OCE,OCB):
  """
  Usage:
    - e.g., Ag,Ac = analytic.rec_lens.qeb(lmax,rlmin,rlmax,fC,OCE,OCB)
  """
  return analytic.rec_lens.qeb(lmax,rlmin,rlmax,fC,OCE,OCB)

def qbb(lmax,rlmin,rlmax,fC,OCB):
  """
  Usage:
    - e.g., Ag,Ac = analytic.rec_lens.qbb(lmax,rlmin,rlmax,fC,OCB)
  """
  return analytic.rec_lens.qbb(lmax,rlmin,rlmax,fC,OCB)

def qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE):
  """
  Usage:
    - e.g., Ig,Ic = analytic.rec_lens.qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE)
  """
  return analytic.rec_lens.qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE)

def qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE):
  """
  Usage:
    - e.g., Ig,Ic = analytic.rec_lens.qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE)
  """
  return analytic.rec_lens.qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE)

def qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE):
  """
  Usage:
    - e.g., Ig,Ic = analytic.rec_lens.qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE)
  """
  return analytic.rec_lens.qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE)

def qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE):
  """
  Usage:
    - e.g., Ig,Ic = analytic.rec_lens.qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE)
  """
  return analytic.rec_lens.qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE)

def qmv(lmax,QDO,Al,Il):
  """
 Compute MV estimator from pre-computed normalization (Al) and correlations (Il)
   set ids
     noise covariance
  Usage:
    - e.g., MV = analytic.rec_lens.qmv(lmax,QDO,Al,Il)
  """
  return analytic.rec_lens.qmv(lmax,QDO,Al,Il)

