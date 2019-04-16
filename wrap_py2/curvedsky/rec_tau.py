import curvedsky

def qtt_sym(lmax,rlmin,rlmax,fC,Tlm,nside=0):
  """
  Reconstructing inhomogeneous tau from the temperature quadratic estimator, assuming Tlm1=Tlm2

  Args:
    - lmax (int)        : maximum multipole of output tau alms
    - rlmin (int)       : minimum multipole of CMB for reconstruction
    - rlmax (int)       : maximum multipole of CMB for reconstruction
    - fC[l] (double)    : temperature angular power spectrum, with bounds (1:rlmax)
    - Tlm[l,m] (dcmplx) : inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)       : Nside for the convolution calculation, default to lmax

  Returns:
    - xlm[l,m] (dcmplx) : inhomogeneous tau alm, with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., xlm = curvedsky.rec_tau.qtt_sym(lmax,rlmin,rlmax,fC,Tlm,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_tau.qtt_sym(lmax,rlmin,rlmax,fC,Tlm,nside)

def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside=0):
  """
  Reconstructing inhomogeneous tau from the temperature quadratic estimator

  Args:
    - lmax (int)        : maximum multipole of output tau alms
    - rlmin (int)       : minimum multipole of CMB for reconstruction
    - rlmax (int)       : maximum multipole of CMB for reconstruction
    - fC[l] (double)    : temperature angular power spectrum, with bounds (1:rlmax)
    - Tlm1[l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    - Tlm2[l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)       : Nside for the convolution calculation, default to lmax

  Returns:
    - xlm[l,m] (dcmplx) : inhomogeneous tau alm, with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., xlm = curvedsky.rec_tau.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_tau.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside)

