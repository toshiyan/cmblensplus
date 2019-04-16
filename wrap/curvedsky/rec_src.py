import curvedsky

def qtt_sym(lmax,rlmin,rlmax,Tlm,nside=0):
  """
  Reconstructing point sources from the temperature quadratic estimator, assuming Tlm1=Tlm2

  Args:
    - lmax (int)        : maximum multipole of output point-source alms
    - rlmin (int)       : minimum multipole of CMB for reconstruction
    - rlmax (int)       : maximum multipole of CMB for reconstruction
    - Tlm[l,m] (dcmplx) : inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)       : Nside for the convolution calculation, default to lmax

  Returns:
    - slm[l,m] (dcmplx) : point-source alm, with bounds (0:lmax,0:lmax)

   alm to map 
   map to alm
  Usage:
    - e.g., slm = curvedsky.rec_src.qtt_sym(lmax,rlmin,rlmax,Tlm,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_src.qtt_sym(lmax,rlmin,rlmax,Tlm,nside)

def qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside=0):
  """
  Reconstructing point sources from the temperature quadratic estimator

  Args:
    - lmax (int)        : maximum multipole of output point-source alms
    - rlmin (int)       : minimum multipole of CMB for reconstruction
    - rlmax (int)       : maximum multipole of CMB for reconstruction
    - Tlm1[l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    - Tlm2[l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)       : Nside for the convolution calculation, default to lmax

  Returns:
    - slm[l,m] (dcmplx) : point-source alm, with bounds (0:lmax,0:lmax)

   alm to map 
   map to alm
  Usage:
    - e.g., slm = curvedsky.rec_src.qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_src.qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside)

