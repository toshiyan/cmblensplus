import libcurvedsky
import numpy

def qeb(lmax,rlmin,rlmax,EB,Elm,Blm,nside_t=0,verbose=False):
  """
  Reconstructing amplitude modulation by the odd EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :EB [*l*] (*double*): EB spectrum, with bounds (0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :verbose (*bool*): Output messages, default to False

  Returns:
    :alm [*l,m*] (*dcmplx*): Reconstructed alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.rec_iamp.qeb(lmax,rlmin,rlmax,EB,Elm,Blm,nside_t,verbose):
  """
  return libcurvedsky.rec_iamp.qeb(lmax,rlmin,rlmax,EB,Elm,Blm,nside_t,verbose)

