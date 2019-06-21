import libflatsky

def qtt(nx,ny,D,rL,OT,TT,eL):
  """
  Normalization of the temperature quadratic estimator for CMB lensing potential and its curl mode

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OT[*lx,ly*] (*double*): Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
    :TT[*lx,ly*] (*double*): Theory temperature power spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :Ag[*lx,ly*] (*dcmplx*): Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
    :Ac[*lx,ly*] (*dcmplx*): Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)

  Usage:
    :Ag,Ac = flatsky.norm_lens.qtt(nx,ny,D,rL,OT,TT,eL):
  """
  return libflatsky.norm_lens.qtt(nx,ny,D,rL,OT,TT,eL)

def qte(nx,ny,D,rL,OT,OE,TE,eL):
  """
  Normalization of the TE quadratic estimator for CMB lensing potential and its curl mode

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OT[*lx,ly*] (*double*): Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
    :OE[*lx,ly*] (*double*): Inverse of Observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
    :TE[*lx,ly*] (*double*): Theory TE cross spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :Ag[*lx,ly*] (*dcmplx*): Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
    :Ac[*lx,ly*] (*dcmplx*): Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)

  Usage:
    :Ag,Ac = flatsky.norm_lens.qte(nx,ny,D,rL,OT,OE,TE,eL):
  """
  return libflatsky.norm_lens.qte(nx,ny,D,rL,OT,OE,TE,eL)

def qtb(nx,ny,D,OT,OB,TE,rL,eL):
  """
  Normalization of the TB quadratic estimator for CMB lensing potential and its curl mode

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OT[*lx,ly*] (*double*): Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
    :OB[*lx,ly*] (*double*): Inverse of Observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
    :TE[*lx,ly*] (*double*): Theory TE cross spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :Ag[*lx,ly*] (*dcmplx*): Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
    :Ac[*lx,ly*] (*dcmplx*): Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)

  Usage:
    :Ag,Ac = flatsky.norm_lens.qtb(nx,ny,D,OT,OB,TE,rL,eL):
  """
  return libflatsky.norm_lens.qtb(nx,ny,D,OT,OB,TE,rL,eL)

def qee(nx,ny,D,OE,EE,rL,eL):
  """
  Normalization of the EE quadratic estimator for CMB lensing potential and its curl mode

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OE[*lx,ly*] (*double*): Inverse of Observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
    :EE[*lx,ly*] (*double*): Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :Ag[*lx,ly*] (*dcmplx*): Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
    :Ac[*lx,ly*] (*dcmplx*): Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)

  Usage:
    :Ag,Ac = flatsky.norm_lens.qee(nx,ny,D,OE,EE,rL,eL):
  """
  return libflatsky.norm_lens.qee(nx,ny,D,OE,EE,rL,eL)

def qeb(nx,ny,D,OE,OB,EE,rL,eL):
  """
  Normalization of the EB quadratic estimator for CMB lensing potential and its curl mode

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :rL[*2*] (*int*): Minimum and maximum multipole of CMB for reconstruction
    :OE[*lx,ly*] (*double*): Inverse of Observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
    :OB[*lx,ly*] (*double*): Inverse of Observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
    :EE[*lx,ly*] (*double*): Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
    :eL[*2*] (*int*): Minimum and maximum multipole of output normalization spectrum, with bounds (2)

  Returns:
    :Ag[*lx,ly*] (*dcmplx*): Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
    :Ac[*lx,ly*] (*dcmplx*): Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)

  Usage:
    :Ag,Ac = flatsky.norm_lens.qeb(nx,ny,D,OE,OB,EE,rL,eL):
  """
  return libflatsky.norm_lens.qeb(nx,ny,D,OE,OB,EE,rL,eL)

