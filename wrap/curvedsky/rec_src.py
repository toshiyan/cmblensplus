import libcurvedsky
import numpy

def qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside_t=0,verbose=False):
  """
  Reconstructing point sources from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output point-source alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :Tlm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Tlm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :verbose (*bool*): Output messages, default to False

  Returns:
    :slm [*l,m*] (*dcmplx*): Point-source alm, with bounds (0:lmax,0:lmax)

  Usage:
    :slm = curvedsky.rec_src.qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside_t,verbose):
  """
  return libcurvedsky.rec_src.qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,nside_t,verbose)

