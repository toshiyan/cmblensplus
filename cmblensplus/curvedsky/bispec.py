import numpy
from ._core import lib_bispec

def make_quad_gauss(alm,bst=1):
  """
  Return a non-Gaussian alm. The corresponding non-Gaussian field is defined as

    delta^NL(n) = delta^L(n) + delta^L(n)**2

  where delta^L(n) is a gaussian field obtained from the input alm.

  Args:
    :alm [l,m] (dcmplx): Input harmonic coefficients, with bounds (0:lmax,0:lmax).

  Returns:
    :qlm [l,m] (dcmplx): Output harmonic coefficients of the non-Gaussian fields, with bounds (0:lmax,0:lmax).

  """
  return lib_bispec.make_quad_gauss(alm,bst=bst)

def bispec_norm(bp,bstype='equi',bst=2,sL=None,nthreads=0):
  """
  Return normalization of the binned reduced bispectrum for a given multipole bin

  Args:
    :bp [edge] (double): Bin edges, with bounds (bn+1)

  Args(optional):
    :bstype (str): Configuration of the bispectrum, default to equilateral
    :bst (int): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    :sL[2] (int): The fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin

  Returns:
    :norm [bin] (double): Normalization of the binned reduced bispectrum at each bin, with bounds (bn)

  """
  return lib_bispec.bispec_norm(bp,bstype=bstype,bst=bst,sL=sL,nthreads=nthreads)

def bispec_bin(bp,alm,bstype='equi',bst=2,sL=None,nthreads=0):
  """
  Return the unnormalized binned reduced bispectrum for a given multipole bin

  Args:
    :bp [edge] (double): Bin edges, with bounds (bn+1)
    :alm [l,m] (dcmplx): Input harmonic coefficients, with bounds (0:lmax,0:lmax)

  Args(optional):
    :bstype (str): Configuration of the bispectrum, default to equilateral
    :bst (int): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    :sL[2] (int): The fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin

  Returns:
    :bis [bin] (double): The unnormalized binned reduced bispectrum at each bin, with bounds (bn)

  """
  return lib_bispec.bispec_bin(bp,alm,bstype,bst,sL,nthreads=nthreads)

def xbispec_bin(bp,alm,bstype='equi',bst=2,sL=None,nthreads=0):
  """
  Return the unnormalized binned reduced cross-bispectrum for a given multipole bin

  Args:
    :bp [edge] (double): Bin edges, with bounds (bn+1)
    :alm [n,l,m] (dcmplx): Input harmonic coefficients, with bounds (0:n-1,0:lmax,0:lmax)

  Args(optional):
    :bstype (str): Configuration of the bispectrum, default to equilateral
    :bst (int): A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
    :sL[2] (int): The fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin

  Returns:
    :bis [bin] (double): The unnormalized binned reduced bispectrum at each bin, with bounds (bn)

  """
  return lib_bispec.xbispec_bin(bp,alm,bstype,bst,sL,nthreads=nthreads)

