import libflatsky

def qeb(nx,ny,D,rL,IE,IB,EE,eL):
  """
  Normalization of the EB quadratic estimator for patchy tau

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :IE[*lx,ly*] (*double*): Inverse of observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
    :IB[*lx,ly*] (*double*): Inverse of observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
    :EE[*lx,ly*] (*double*): Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :At[*lx,ly*] (*dcmplx*): Normalization of patchy tau on 2D grid, with bounds (nx,ny)

  Usage:
    :At = flatsky.norm_tau.qeb(nx,ny,D,rL,IE,IB,EE,eL):
  """
  return libflatsky.norm_tau.qeb(nx,ny,D,rL,IE,IB,EE,eL)

