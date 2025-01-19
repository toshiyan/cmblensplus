from cmblensplus import libcurvedsky
import numpy

def make_quad_gauss(alm):
  """
  Return a non-Gaussian alm. The corresponding non-Gaussian field is defined as

    delta^NL(n) = delta^L(n) + delta^L(n)**2

  where delta^L(n) is a gaussian field obtained from the input alm.

  Args:
    :alm [*l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Returns:
    :qlm [*l,m*] (*dcmplx*): Output harmonic coefficients of the non-Gaussian fields, with bounds (0:lmax,0:lmax).

  Usage:
    :qlm = curvedsky.bispec.make_quad_gauss(lmax,alm):
  """
  lmax = len(alm[:,0]) - 1
  return libcurvedsky.bispec.make_quad_gauss(lmax,alm)

def bispec_norm(bp,bstype='equi',bst=2,sL=None):
  """
  Return normalization of the binned reduced bispectrum for a given multipole bin

  Args:
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
  if sL is None: sL = numpy.zeros(2)
  bn = len(bp) - 1
  return libcurvedsky.bispec.bispec_norm(bn,bp,bstype,bst,sL)

def bispec_bin(bp,alm,bstype='equi',bst=2,sL=None):
  """
  Return the unnormalized binned reduced bispectrum for a given multipole bin

  Args:
    :bp [*edge*] (*double*): Bin edges, with bounds (bn+1)
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
  if sL is None: sL = numpy.zeros(2)
  bn = len(bp) - 1
  lmax = len(alm[:,0]) - 1
  return libcurvedsky.bispec.bispec_bin(bn,bp,lmax,alm,bstype,bst,sL)

def xbispec_bin(bp,alm,bstype='equi',bst=2,sL=None):
  """
  Return the unnormalized binned reduced cross-bispectrum for a given multipole bin

  Args:
    :bp [*edge*] (*double*): Bin edges, with bounds (bn+1)
    :alm [*n,l,m*] (*dcmplx*): Input harmonic coefficients, with bounds (0:n-1,0:lmax,0:lmax)

  Args(optional):
    :bstype (*str*): Configuration of the bispectrum, default to equilateral
    :bst (*int*): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    :sL[*2*] (*int*): The fixed bin for the squeezed configuration, b[*sL,eL,eL*], default to the lowest multipole bin

  Returns:
    :bis [*bin*] (*double*): The unnormalized binned reduced bispectrum at each bin, with bounds (bn)

  Usage:
    :bis = curvedsky.bispec.xbispec_bin(bn,bp,lmax,n,alm,bstype,bst,sL):
  """
  if sL is None: sL = numpy.zeros(2)
  bn = len(bp) - 1
  n = len(alm[:,0,0])
  lmax = len(alm[0,:,0]) - 1
  return libcurvedsky.bispec.xbispec_bin(bn,bp,lmax,n,alm,bstype,bst,sL)

