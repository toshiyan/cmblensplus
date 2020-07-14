import libflatsky

def qtt(nx,ny,D,rL,OT,TT,eL):
  """
  Normalization of the temperature quadratic estimators between CMB lensing potential and patchy tau

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OT[*lx,ly*] (*double*): Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
    :TT[*lx,ly*] (*double*): Theory temperature power spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :Aks[*lx,ly*] (*dcmplx*): Lensing-tau cross normalization on 2D grid, with bounds (nx,ny)

  Usage:
    :Aks = flatsky.norm_kxs.qtt(nx,ny,D,rL,OT,TT,eL):
  """
  return libflatsky.norm_kxs.qtt(nx,ny,D,rL,OT,TT,eL)

