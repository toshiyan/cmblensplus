from cmblensplus import libflatsky
import numpy

def qeb(nx,ny,D,rL,IE,IB,EE,eL,BB=0):
  """
  Normalization of the EB quadratic estimator for anisotropic pol. rot. angles

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :IE[*lx,ly*] (*double*): Inverse of observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
    :IB[*lx,ly*] (*double*): Inverse of observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
    :EE[*lx,ly*] (*double*): Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Args(Optional):
    :BB[*lx,ly*] (*double*): Theory B-mode spectrum on 2D grid, with bounds (nx,ny), default to BB=0

  Returns:
    :Aa[*lx,ly*] (*dcmplx*): Normalization of anisotropic pol. rot. angles on 2D grid, with bounds (nx,ny)

  Usage:
    :Aa = flatsky.norm_rot.qeb(nx,ny,D,rL,IE,IB,EE,eL,BB):
  """
  return libflatsky.norm_rot.qeb(nx,ny,D,rL,IE,IB,EE,eL,BB)

