import libcurvedsky

def qtt_sym(lmax,rlmin,rlmax,fC,Tlm,nside=None):
  """
  Reconstructing inhomogeneous tau from the temperature quadratic estimator, assuming Tlm1=Tlm2

  Args:
    :lmax (*int*): Maximum multipole of output tau alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Temperature angular power spectrum, with bounds (1:rlmax)
    :Tlm [*l,m*] (*dcmplx*): Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside (*int*): Nside for the convolution calculation, default to lmax

  Returns:
    :xlm [*l,m*] (*dcmplx*): Inhomogeneous tau alm, with bounds (0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.rec_tau.qtt_sym(lmax,rlmin,rlmax,fC,Tlm,nside):
  """
  if nside is None: nside= lmax
  return libcurvedsky.rec_tau.qtt_sym(lmax,rlmin,rlmax,fC,Tlm,nside)

def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside=None):
  """
  Reconstructing inhomogeneous tau from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output tau alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): Temperature angular power spectrum, with bounds (1:rlmax)
    :Tlm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Tlm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside (*int*): Nside for the convolution calculation, default to lmax

  Returns:
    :xlm [*l,m*] (*dcmplx*): Inhomogeneous tau alm, with bounds (0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.rec_tau.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside):
  """
  if nside is None: nside= lmax
  return libcurvedsky.rec_tau.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside)

def oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,nside=None):
  """
  Reconstructing amplitude modulation by the odd EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fEB [*l*] (*double*): EB spectrum, with bounds (0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside (*int*): Nside for the convolution calculation, default to lmax

  Returns:
    :alm [*l,m*] (*dcmplx*): Reconstructed alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.rec_tau.oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,nside):
  """
  if nside is None: nside= lmax
  return libcurvedsky.rec_tau.oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,nside)

