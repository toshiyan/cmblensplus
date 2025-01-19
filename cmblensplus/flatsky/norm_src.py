from cmblensplus import libflatsky
import numpy

def qtt(nx,ny,D,rL,OT,eL):
  """
  Normalization of the temperature quadratic estimator for point source

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OT[*lx,ly*] (*double*): Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :As[*lx,ly*] (*dcmplx*): Normalization of point source on 2D grid, with bounds (nx,ny)

  Usage:
    :As = flatsky.norm_src.qtt(nx,ny,D,rL,OT,eL):
  """
  return libflatsky.norm_src.qtt(nx,ny,D,rL,OT,eL)

def qeb(nx,ny,D,rL,IE,IB,eL):
  """
  Normalization of the EB quadratic estimator for point source 

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :IE[*lx,ly*] (*double*): Inverse of observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
    :IB[*lx,ly*] (*double*): Inverse of observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :As[*lx,ly*] (*dcmplx*): Normalization of point source on 2D grid, with bounds (nx,ny)

  Usage:
    :As = flatsky.norm_src.qeb(nx,ny,D,rL,IE,IB,eL):
  """
  return libflatsky.norm_src.qeb(nx,ny,D,rL,IE,IB,eL)

