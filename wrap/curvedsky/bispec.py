import libcurvedsky

def make_quad_gauss(lmax,alm):
  """
  Return a non-Gaussian alm. The corresponding non-Gaussian field is defined as

    delta^NL(n) = delta^L(n) + delta^L(n)**2

  where delta^L(n) is a gaussian field obtained from the input alm.

  Args:
    :lmax (*int*): Maximum multipole of alm
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Returns:
    :qlm [*l,m*] (*dcmplx*): Output harmonic coefficients of the non-Gaussian fields, with bounds (0:lmax,0:lmax).

  Usage:
    :qlm = curvedsky.bispec.make_quad_gauss(lmax,alm):
  """
  return libcurvedsky.bispec.make_quad_gauss(lmax,alm)

def bispec_norm(bn,bp,bstype= 'equi',bst= 2,sL=None):
  """
  Return normalization of the binned reduced bispectrum for a given multipole bin

  Args:
    :bn (*int*): Number of multipole bins
    :bp [*edge*] (*double*): Bin edges, with bounds (bn+1)

  Args(optional):
    :bstype (*str*): Configuration of the bispectrum, default to equilateral
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    :sL[*2*] (*int*): The fixed bin for the squeezed configuration, b[*sL,eL,eL*], default to the lowest multipole bin

  Returns:
    :norm [*bin*] (*double*): Normalization of the binned reduced bispectrum at each bin, with bounds (bn)

  Usage:
    :norm = curvedsky.bispec.bispec_norm(bn,bp,bstype,bst,sL):
  """
  if sL is None: sL= [int(bp[0]),int(bp[1])]
  return libcurvedsky.bispec.bispec_norm(bn,bp,bstype,bst,sL)

def bispec_bin(bn,bp,lmax,alm,bstype= 'equi',bst= 'equi',sL=None):
  """
  Return the unnormalized binned reduced bispectrum for a given multipole bin

  Args:
    :bn (*int*): Number of multipole bins
    :bp [*edge*] (*double*): Bin edges, with bounds (bn+1)
    :lmax (*int*): Maximum multipole of the input alm
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax)

  Args(optional):
    :bstype (*str*): Configuration of the bispectrum, default to equilateral
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    :sL[*2*] (*int*): The fixed bin for the squeezed configuration, b[*sL,eL,eL*], default to the lowest multipole bin

  Returns:
    :bis [*bin*] (*double*): The unnormalized binned reduced bispectrum at each bin, with bounds (bn)

  Usage:
    :bis = curvedsky.bispec.bispec_bin(bn,bp,lmax,alm,bstype,bst,sL):
  """
  if sL is None: sL= [int(bp[0]),int(bp[1])]
  return libcurvedsky.bispec.bispec_bin(bn,bp,lmax,alm,bstype,bst,sL)

def equi(lmin,lmax,alm,bst= 2):
  """
  Compute equilateral shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l,l]

  Args:
    :lmin (*int*): Minimum multipole of the bin
    :lmax (*int*): Maximum multipole of the bin
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    :bispec (*double*): Unnormalized binned reduced bispectrum at the bin, [*lmin,lmax*]

  Usage:
    :bispec = curvedsky.bispec.equi(lmin,lmax,alm,bst):
  """
  return libcurvedsky.bispec.equi(lmin,lmax,alm,bst)

def fold(lmin,lmax,alm,bst= 2):
  """
  Compute folded shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l/2,l/2]

  Args:
    :lmin (*int*): Minimum multipole of the bin
    :lmax (*int*): Maximum multipole of the bin
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    :bispec (*double*): Unnormalized binned reduced bispectrum at the bin, [*lmin,lmax*]

  Usage:
    :bispec = curvedsky.bispec.fold(lmin,lmax,alm,bst):
  """
  return libcurvedsky.bispec.fold(lmin,lmax,alm,bst)

def sque(eL,sL,l1,alm,bst= 2):
  """
  Compute squeezed shape of the unnormalized binned reduced bispectrum for a given alm, b[sL,eL,eL]

  Args:
    :eL[*2*] (*int*): Minimum and maximum multipoles of the bin, with bounds (2)
    :sL[*2*] (*int*): Minimum and maximum multipoles of the fixed bin, with bounds (2)
    :l1 (*int*): Maximum multipole of the input alm, satisfying eLmax,sLmax<=l1
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    :bispec (*double*): Unnormalized binned reduced bispectrum at the bin, [*lmin,lmax*]

  Usage:
    :bispec = curvedsky.bispec.sque(eL,sL,l1,alm,bst):
  """
  return libcurvedsky.bispec.sque(eL,sL,l1,alm,bst)

def isos(eL,aL,l1,alm,bst= 2):
  """
  Compute isosceles shape of the unnormalized binned reduced bispectrum for a given alm, b[eL,aL,aL]

  Args:
    :eL[*2*] (*int*): Minimum and maximum multipoles of the bin, with bounds (2)
    :aL[*2*] (*int*): Minimum and maximum multipoles of the fixed bin, with bounds (2)
    :l1 (*int*): Maximum multipole of the input alm, satisfying eLmax,sLmax<=l1
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    :bispec (*double*): Unnormalized binned reduced bispectrum at the bin, [*lmin,lmax*]

  Usage:
    :bispec = curvedsky.bispec.isos(eL,aL,l1,alm,bst):
  """
  return libcurvedsky.bispec.isos(eL,aL,l1,alm,bst)

