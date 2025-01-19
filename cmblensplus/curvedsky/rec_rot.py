from cmblensplus import libcurvedsky
import numpy

def qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,nside_t=0,verbose=False):
  """
  Reconstructing pol rotation angle from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCE [*l*] (*double*): EE spectrum, with bounds (0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :verbose (*bool*): Output messages, default to False

  Returns:
    :alm [*l,m*] (*dcmplx*): Rotation angle alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.rec_rot.qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,nside_t,verbose):
  """
  return libcurvedsky.rec_rot.qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,nside_t,verbose)

