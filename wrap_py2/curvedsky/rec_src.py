import libcurvedsky

def qtt_sym(lmax,rlmin,rlmax,Tlm,nside=None):
  """
  Reconstructing point sources from the temperature quadratic estimator, assuming Tlm1=Tlm2

  Args:
    :lmax (*int*): Maximum multipole of output point-source alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :Tlm [*l,m*] (*dcmplx*): Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside (*int*): Nside for the convolution calculation, default to lmax

  Returns:
    :slm [*l,m*] (*dcmplx*): Point-source alm, with bounds (0:lmax,0:lmax)

  Usage:
    :slm = curvedsky.rec_src.qtt_sym(lmax,rlmin,rlmax,Tlm,nside):
  """
  if nside is None: nside= lmax
  return libcurvedsky.rec_src.qtt_sym(lmax,rlmin,rlmax,Tlm,nside)

def qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside=None):
  """
  Reconstructing point sources from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output point-source alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :Tlm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Tlm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside (*int*): Nside for the convolution calculation, default to lmax

  Returns:
    :slm [*l,m*] (*dcmplx*): Point-source alm, with bounds (0:lmax,0:lmax)

  Usage:
    :slm = curvedsky.rec_src.qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside):
  """
  if nside is None: nside= lmax
  return libcurvedsky.rec_src.qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside)

