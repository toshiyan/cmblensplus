import libflatsky

def qtt(nx,ny,D,rL,fC,T1,T2):
  """
  Reconstructing patchy tau from the temperature quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): Temperature power spectrum on 2D grid, with bounds (nx,ny)
    :T1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered temperature, with bounds (nx,ny)
    :T2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered temperature, with bounds (nx,ny)

  Returns:
    :tlm[*lx,ly*] (*dcmplx*): 2D Fourier modes of patchy tau, with bounds (nx,ny)

  Usage:
    :tlm = flatsky.rec_tau.qtt(nx,ny,D,rL,fC,T1,T2):
  """
  return libflatsky.rec_tau.qtt(nx,ny,D,rL,fC,T1,T2)

def qte(nx,ny,D,rL,fC,T,E):
  """
  Reconstructing patchy tau from the suboptimal TE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)

  Returns:
    :tlm[*lx,ly*] (*dcmplx*): 2D Fourier modes of patchy tau, with bounds (nx,ny)

  Usage:
    :tlm = flatsky.rec_tau.qte(nx,ny,D,rL,fC,T,E):
  """
  return libflatsky.rec_tau.qte(nx,ny,D,rL,fC,T,E)

def qtb(nx,ny,D,rL,fC,T,B):
  """
  Reconstructing patchy tau from the TB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): TE cross power spectrum on 2D grid, with bounds (nx,ny)
    :T[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Returns:
    :tlm[*lx,ly*] (*dcmplx*): 2D Fourier modes of patchy tau, with bounds (nx,ny)

  Usage:
    :tlm = flatsky.rec_tau.qtb(nx,ny,D,rL,fC,T,B):
  """
  return libflatsky.rec_tau.qtb(nx,ny,D,rL,fC,T,B)

def qee(nx,ny,D,rL,fC,E1,E2):
  """
  Reconstructing patchy tau from the EE quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fC[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :E1[*lx,ly*] (*dcmplx*): 2D Fourier modes of 1st inverse-variance filtered E-mode, with bounds (nx,ny)
    :E2[*lx,ly*] (*dcmplx*): 2D Fourier modes of 2nd inverse-variance filtered E-mode, with bounds (nx,ny)

  Returns:
    :tlm[*lx,ly*] (*dcmplx*): 2D Fourier modes of patchy tau, with bounds (nx,ny)

  Usage:
    :tlm = flatsky.rec_tau.qee(nx,ny,D,rL,fC,E1,E2):
  """
  return libflatsky.rec_tau.qee(nx,ny,D,rL,fC,E1,E2)

def qeb(nx,ny,D,rL,fE,fB,E,B):
  """
  Reconstructing patchy tau from the EB quadratic estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :fE[*lx,ly*] (*double*): EE power spectrum on 2D grid, with bounds (nx,ny)
    :fB[*lx,ly*] (*double*): BB power spectrum on 2D grid, with bounds (nx,ny)
    :E[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
    :B[*lx,ly*] (*dcmplx*): 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)

  Returns:
    :tlm[*lx,ly*] (*dcmplx*): 2D Fourier modes of patchy tau, with bounds (nx,ny)

  Usage:
    :tlm = flatsky.rec_tau.qeb(nx,ny,D,rL,fE,fB,E,B):
  """
  return libflatsky.rec_tau.qeb(nx,ny,D,rL,fE,fB,E,B)

