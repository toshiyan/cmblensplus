import libflatsky
import numpy

def qtt(nx,ny,D,rL,T1,T2):
  """
  Reconstructing point source fields from the temperature quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): Temperature power spectrum on 2D grid, with bounds (nx,ny)
    :T1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered temperature, with bounds (nx,ny)
    :T2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered temperature, with bounds (nx,ny)

  Returns:
    :slm[*lx,ly*] (*dcmplx*): 2D Fourier modes of point source fields, with bounds (nx,ny)

  Usage:
    :slm = flatsky.rec_src.qtt(nx,ny,D,rL,T1,T2):
  """
  return libflatsky.rec_src.qtt(nx,ny,D,rL,T1,T2)

def qte(nx,ny,D,rL,T,E):
  """
  Reconstructing point source fields from the suboptimal TE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)

  Returns:
    :slm[*lx,ly*] (*dcmplx*): 2D Fourier modes of point source fields, with bounds (nx,ny)

  Usage:
    :slm = flatsky.rec_src.qte(nx,ny,D,rL,T,E):
  """
  return libflatsky.rec_src.qte(nx,ny,D,rL,T,E)

def qtb(nx,ny,D,rL,T,B):
  """
  Reconstructing point source fields from the TB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Returns:
    :slm[*lx,ly*] (*dcmplx*): 2D Fourier modes of point source fields, with bounds (nx,ny)

  Usage:
    :slm = flatsky.rec_src.qtb(nx,ny,D,rL,T,B):
  """
  return libflatsky.rec_src.qtb(nx,ny,D,rL,T,B)

def qee(nx,ny,D,rL,E1,E2):
  """
  Reconstructing point source fields from the EE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered E-mode, with bounds (nx,ny)
    :E2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered E-mode, with bounds (nx,ny)

  Returns:
    :slm[*lx,ly*] (*dcmplx*): 2D Fourier modes of point source fields, with bounds (nx,ny)

  Usage:
    :slm = flatsky.rec_src.qee(nx,ny,D,rL,E1,E2):
  """
  return libflatsky.rec_src.qee(nx,ny,D,rL,E1,E2)

def qeb(nx,ny,D,rL,E,B):
  """
  Reconstructing point source fields from the EB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Returns:
    :slm[*lx,ly*] (*dcmplx*): 2D Fourier modes of point source fields, with bounds (nx,ny)

  Usage:
    :slm = flatsky.rec_src.qeb(nx,ny,D,rL,E,B):
  """
  return libflatsky.rec_src.qeb(nx,ny,D,rL,E,B)

