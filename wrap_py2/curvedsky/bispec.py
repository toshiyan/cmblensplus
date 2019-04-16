import curvedsky

def make_quad_gauss(lmax,alm):
  """
  Return a non-Gaussian alm. The corresponding non-Gaussian field is defined as

    delta^NL(n) = delta^L(n) + delta^L(n)**2

  where delta^L(n) is a gaussian field obtained from the input alm.

  Args:
    - lmax (int)        : maximum multipole of alm
    - alm[l,m] (dcmplx) : input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Returns:
    - qlm[l,m] (dcmplx) : output harmonic coefficients of the non-Gaussian fields, with bounds (0:lmax,0:lmax).

  Usage:
    - e.g., qlm = curvedsky.bispec.make_quad_gauss(lmax,alm)
  """
  return curvedsky.bispec.make_quad_gauss(lmax,alm)

def bispec_norm(bn,bp,bstype=0,bst=0,sL=0):
  """
  Return normalization of the binned reduced bispectrum for a given multipole bin

  Args:
    - bn (int)          : number of multipole bins
    - bp[edge] (double) : bin edges, with bounds (bn+1)

  Args(optional):
    - bstype (str) : configuration of the bispectrum, default to equilateral
    - bst (int)    : a parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    - sL[2] (int)  : the fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin

  Returns:
    - norm[bin] (double) : normalization of the binned reduced bispectrum at each bin, with bounds (bn)

  Usage:
    - e.g., norm = curvedsky.bispec.bispec_norm(bn,bp,bstype,bst,sL)
  """
  if bstype==0: bstype= 'equi'
  if bst==0: bst= 2
  if sL==0: sL= 0
  return curvedsky.bispec.bispec_norm(bn,bp,bstype,bst,sL)

def bispec_bin(bn,bp,lmax,alm,bstype=0,bst=0,sL=0):
  """
  Return the unnormalized binned reduced bispectrum for a given multipole bin

  Args:
    - bn (int)          : number of multipole bins
    - bp[edge] (double) : bin edges, with bounds (bn+1)
    - lmax (int)        : maximum multipole of the input alm
    - alm[l,m] (dcmplx) : input harmonic coefficients, with bounds (0:lmax,0:lmax)

  Args(optional):
    - bstype (str) : configuration of the bispectrum, default to equilateral
    - bst (int) : a parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    - sL[2] (int)  : the fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin

  Returns:
    - bis[bin] (double) : the unnormalized binned reduced bispectrum at each bin, with bounds (bn)

  Usage:
    - e.g., bis = curvedsky.bispec.bispec_bin(bn,bp,lmax,alm,bstype,bst,sL)
  """
  if bstype==0: bstype= 'equi'
  if bst==0: bst= 'equi'
  if sL==0: sL= 0
  return curvedsky.bispec.bispec_bin(bn,bp,lmax,alm,bstype,bst,sL)

def equi(lmin,lmax,alm,bst=0):
  """
  Compute equilateral shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l,l]

  Args:
    - lmin (int)        : minimum multipole of the bin
    - lmax (int)        : maximum multipole of the bin
    - alm[l,m] (dcmplx) : input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    - bst (int)         : a parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    - bispec (double)   : unnormalized binned reduced bispectrum at the bin, [lmin,lmax]

  Usage:
    - e.g., bispec = curvedsky.bispec.equi(lmin,lmax,alm,bst)
  """
  if bst==0: bst= 2
  return curvedsky.bispec.equi(lmin,lmax,alm,bst)

def fold(lmin,lmax,alm,bst=0):
  """
  Compute folded shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l/2,l/2]

  Args:
    - lmin (int)        : minimum multipole of the bin
    - lmax (int)        : maximum multipole of the bin
    - alm[l,m] (dcmplx) : input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    - bst (int)         : a parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    - bispec (double)   : unnormalized binned reduced bispectrum at the bin, [lmin,lmax]

  Usage:
    - e.g., bispec = curvedsky.bispec.fold(lmin,lmax,alm,bst)
  """
  if bst==0: bst= 2
  return curvedsky.bispec.fold(lmin,lmax,alm,bst)

def sque(eL,sL,l1,alm,bst=0):
  """
  Compute squeezed shape of the unnormalized binned reduced bispectrum for a given alm, b[sL,eL,eL]

  Args:
    - eL[] (int)        : minimum and maximum multipoles of the bin, with bounds (2)
    - sL[] (int)        : minimum and maximum multipoles of the fixed bin, with bounds (2)
    - l1 (int)          : maximum multipole of the input alm, satisfying eLmax,sLmax<=l1
    - alm[l,m] (dcmplx) : input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    - bst (int)         : a parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    - bispec (double)   : unnormalized binned reduced bispectrum at the bin, [lmin,lmax]

  Usage:
    - e.g., bispec = curvedsky.bispec.sque(eL,sL,l1,alm,bst)
  """
  if bst==0: bst= 2
  return curvedsky.bispec.sque(eL,sL,l1,alm,bst)

def isos(eL,aL,l1,alm,bst=0):
  """
  Compute isosceles shape of the unnormalized binned reduced bispectrum for a given alm, b[eL,aL,aL]

  Args:
    - eL[] (int)        : minimum and maximum multipoles of the bin, with bounds (2)
    - aL[] (int)        : minimum and maximum multipoles of the fixed bin, with bounds (2)
    - l1 (int)          : maximum multipole of the input alm, satisfying eLmax,sLmax<=l1
    - alm[l,m] (dcmplx) : input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Args(optional):
    - bst (int)         : a parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.

  Returns:
    - bispec (double)   : unnormalized binned reduced bispectrum at the bin, [lmin,lmax]

  Usage:
    - e.g., bispec = curvedsky.bispec.isos(eL,aL,l1,alm,bst)
  """
  if bst==0: bst= 2
  return curvedsky.bispec.isos(eL,aL,l1,alm,bst)

