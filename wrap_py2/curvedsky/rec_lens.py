import curvedsky

def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside=0):
  """
  Reconstructing CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    - lmax (int)        : maximum multipole of output lensing potential alms
    - rlmin (int)       : minimum multipole of CMB for reconstruction
    - rlmax (int)       : maximum multipole of CMB for reconstruction
    - fC[l] (double)    : temperature angular power spectrum, with bounds (0:rlmax)
    - Tlm1[l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    - Tlm2[l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)       : Nside for the convolution calculation, default to lmax

  Returns:
    - glm[l,m] (dcmplx) : CMB lensing potential alm, with bounds (0:lmax,0:lmax)
    - clm[l,m] (dcmplx) : curl mode (pseudo lensing potential) alm, with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., glm,clm = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside)

def qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside=0):
  """
  Reconstructing CMB lensing potential and its curl mode from the TE quadratic estimator

  Args:
    - lmax (int)   : maximum multipole of output lensing potential alms
    - rlmin (int)  : minimum multipole of CMB for reconstruction
    - rlmax (int)  : maximum multipole of CMB for reconstruction
    - fC[l] (double): temperature-E cross-angular power spectrum, with bounds (0:rlmax)
    - Tlm  (dcmplx): inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    - Elm  (dcmplx): inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)  : Nside for the convolution calculation, default to lmax

  Returns:
    - glm (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
    - clm (dcmplx) : curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., glm,clm = curvedsky.rec_lens.qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_lens.qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside)

def qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside=0):
  """
  Reconstructing CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    - lmax (int)   : maximum multipole of output lensing potential alms
    - rL   (int)   : minimum rL[0] and maximum multipoles rlmax of CMB, with bounds (2)
    - fC[l] (double): temperature-E cross-angular power spectrum, with bounds (0:rlmax)
    - Tlm  (dcmplx): inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    - Blm  (dcmplx): inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)  : Nside for the convolution calculation, default to lmax

  Returns:
    - glm (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
    - clm (dcmplx) : curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., glm,clm = curvedsky.rec_lens.qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_lens.qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside)

def qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside=0):
  """
  Reconstructing CMB lensing potential and its curl mode from the EE quadratic estimator

  Args:
    - lmax (int)   : maximum multipole of output lensing potential alms
    - rL   (int)   : minimum rL[0] and maximum multipoles rlmax of CMB, with bounds (2)
    - fC[l] (double): E-mode angular power spectrum, with bounds (0:rlmax)
    - Elm1 (dcmplx): 1st inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    - Elm2 (dcmplx): 2nd inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)  : Nside for the convolution calculation, default to lmax

  Returns:
    - glm (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
    - clm (dcmplx) : curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., glm,clm = curvedsky.rec_lens.qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_lens.qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside)

def qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside=0):
  """
  Reconstructing CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    - lmax (int)   : maximum multipole of output lensing potential alms
    - rL   (int)   : minimum rL[0] and maximum multipoles rlmax of CMB, with bounds (2)
    - fC[l] (double): E-mode angular power spectrum, with bounds (0:rlmax)
    - Elm  (dcmplx): inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    - Blm  (dcmplx): inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)  : Nside for the convolution calculation, default to lmax

  Returns:
    - glm (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
    - clm (dcmplx) : curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., glm,clm = curvedsky.rec_lens.qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_lens.qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside)

def qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside=0):
  """
  Reconstructing CMB lensing potential and its curl mode from the BB quadratic estimator

  Args:
    - lmax (int)   : maximum multipole of output lensing potential alms
    - rL   (int)   : minimum rL[0] and maximum multipoles rlmax of CMB, with bounds (2)
    - fC[l] (double): B-mode angular power spectrum, with bounds (0:rlmax)
    - Blm1 (dcmplx): 1st inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
    - Blm2 (dcmplx): 2nd inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    - nside (int)  : Nside for the convolution calculation, default to lmax

  Returns:
    - glm (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
    - clm (dcmplx) : curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., glm,clm = curvedsky.rec_lens.qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.rec_lens.qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside)

