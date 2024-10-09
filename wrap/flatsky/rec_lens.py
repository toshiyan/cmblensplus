import libflatsky
import numpy

def qtt(nx,ny,D,rL,fC,T1,T2,gtype= ''):
  """
  Reconstructing CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): Temperature power spectrum on 2D grid, with bounds (nx,ny)
    :T1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered temperature, with bounds (nx,ny)
    :T2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered temperature, with bounds (nx,ny)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :glm[*lx,ly*] (*dcmplx*): 2D Fourier modes of CMB lensing potential, with bounds (nx,ny)
    :clm[*lx,ly*] (*dcmplx*): 2D Fourier modes of Curl mode (pseudo lensing potential), with bounds (nx,ny)

  Usage:
    :glm,clm = flatsky.rec_lens.qtt(nx,ny,D,rL,fC,T1,T2,gtype):
  """
  return libflatsky.rec_lens.qtt(nx,ny,D,rL,fC,T1,T2,gtype)

def qte(nx,ny,D,rL,fC,T,E,gtype= ''):
  """
  Reconstructing CMB lensing potential and its curl mode from the suboptimal TE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :glm[*lx,ly*] (*dcmplx*): 2D Fourier modes of CMB lensing potential, with bounds (nx,ny)
    :clm[*lx,ly*] (*dcmplx*): 2D Fourier modes of Curl mode (pseudo lensing potential), with bounds (nx,ny)

  Usage:
    :glm,clm = flatsky.rec_lens.qte(nx,ny,D,rL,fC,T,E,gtype):
  """
  return libflatsky.rec_lens.qte(nx,ny,D,rL,fC,T,E,gtype)

def qtb(nx,ny,D,rL,fC,T,B,gtype= ''):
  """
  Reconstructing CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :glm[*lx,ly*] (*dcmplx*): 2D Fourier modes of CMB lensing potential, with bounds (nx,ny)
    :clm[*lx,ly*] (*dcmplx*): 2D Fourier modes of Curl mode (pseudo lensing potential), with bounds (nx,ny)

  Usage:
    :glm,clm = flatsky.rec_lens.qtb(nx,ny,D,rL,fC,T,B,gtype):
  """
  return libflatsky.rec_lens.qtb(nx,ny,D,rL,fC,T,B,gtype)

def qee(nx,ny,D,rL,fC,E1,E2,gtype= ''):
  """
  Reconstructing CMB lensing potential and its curl mode from the EE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered E-mode, with bounds (nx,ny)
    :E2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered E-mode, with bounds (nx,ny)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :glm[*lx,ly*] (*dcmplx*): 2D Fourier modes of CMB lensing potential, with bounds (nx,ny)
    :clm[*lx,ly*] (*dcmplx*): 2D Fourier modes of Curl mode (pseudo lensing potential), with bounds (nx,ny)

  Usage:
    :glm,clm = flatsky.rec_lens.qee(nx,ny,D,rL,fC,E1,E2,gtype):
  """
  return libflatsky.rec_lens.qee(nx,ny,D,rL,fC,E1,E2,gtype)

def qeb(nx,ny,D,rL,fC,E,B,gtype= ''):
  """
  Reconstructing CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :glm[*lx,ly*] (*dcmplx*): 2D Fourier modes of CMB lensing potential, with bounds (nx,ny)
    :clm[*lx,ly*] (*dcmplx*): 2D Fourier modes of Curl mode (pseudo lensing potential), with bounds (nx,ny)

  Usage:
    :glm,clm = flatsky.rec_lens.qeb(nx,ny,D,rL,fC,E,B,gtype):
  """
  return libflatsky.rec_lens.qeb(nx,ny,D,rL,fC,E,B,gtype)

def qbb(nx,ny,D,rL,fC,B1,B2,gtype= ''):
  """
  Reconstructing CMB lensing potential and its curl mode from the BB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): BB power spectrum on 2D grid, with bounds (nx,ny)
    :B1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered B-mode, with bounds (nx,ny)
    :B2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered B-mode, with bounds (nx,ny)

  Args(optional):
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)

  Returns:
    :glm[*lx,ly*] (*dcmplx*): 2D Fourier modes of CMB lensing potential, with bounds (nx,ny)
    :clm[*lx,ly*] (*dcmplx*): 2D Fourier modes of Curl mode (pseudo lensing potential), with bounds (nx,ny)

  Usage:
    :glm,clm = flatsky.rec_lens.qbb(nx,ny,D,rL,fC,B1,B2,gtype):
  """
  return libflatsky.rec_lens.qbb(nx,ny,D,rL,fC,B1,B2,gtype)

