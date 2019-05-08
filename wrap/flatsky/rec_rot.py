import libflatsky

def qte(nx,ny,D,rL,fC,T,E):
  """
  Reconstructing anisotropic pol. rot. angles from TE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)

  Returns:
    :alm[*lx,ly*] (*dcmplx*): 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)

  Usage:
    :alm = flatsky.rec_rot.qte(nx,ny,D,rL,fC,T,E):
  """
  return libflatsky.rec_rot.qte(nx,ny,D,rL,fC,T,E)

def qtb(nx,ny,D,rL,fC,T,B):
  """
  Reconstructing anisotropic pol. rot. angles from TB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Returns:
    :alm[*lx,ly*] (*dcmplx*): 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)

  Usage:
    :alm = flatsky.rec_rot.qtb(nx,ny,D,rL,fC,T,B):
  """
  return libflatsky.rec_rot.qtb(nx,ny,D,rL,fC,T,B)

def qee(nx,ny,D,rL,fC,E1,E2):
  """
  Reconstructing anisotropic pol. rot. angles from EE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered E-mode, with bounds (nx,ny)
    :E2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered E-mode, with bounds (nx,ny)

  Returns:
    :alm[*lx,ly*] (*dcmplx*): 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)

  Usage:
    :alm = flatsky.rec_rot.qee(nx,ny,D,rL,fC,E1,E2):
  """
  return libflatsky.rec_rot.qee(nx,ny,D,rL,fC,E1,E2)

def qeb(nx,ny,D,rL,EE,E,B,BB= 0):
  """
  Reconstructing anisotropic pol. rot. angles from EB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Args(Optional):
    :BB[*lx,ly*] (*double*): Theory B-mode spectrum on 2D grid, with bounds (nx,ny), default to BB=0

  Returns:
    :alm[*lx,ly*] (*dcmplx*): 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)

  Usage:
    :alm = flatsky.rec_rot.qeb(nx,ny,D,rL,EE,E,B,BB):
  """
  return libflatsky.rec_rot.qeb(nx,ny,D,rL,EE,E,B,BB)

