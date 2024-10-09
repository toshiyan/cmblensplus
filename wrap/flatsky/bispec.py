import libflatsky
import numpy

def bispec_norm(nx,ny,D,bp,bn=1,dbin_max=-1):
  """
  Usage:
    :norm = flatsky.bispec.bispec_norm(nx,ny,D,bp,dbin_max,bn):
  """
  bn = len(bp) - 1
  if dbmax==-1: dbin_max = bn
  return libflatsky.bispec.bispec_norm(nx,ny,D,bp,dbin_max,bn)

def bispec_bin(kmap,bp,bn=1,kn=1,nx=0,ny=0,dbin_max=-1):
  """
  Usage:
    :bispec = flatsky.bispec.bispec_bin(kn,bn,nx,ny,kmap,bp,dbin_max):
  """
  bn = len(bp) - 1
  kn = len(kmap[:,0,0,0])
  nx = len(kmap[0,0,:,0])
  ny = len(kmap[0,0,0,:])
  if dbmax==-1: dbin_max = bn
  return libflatsky.bispec.bispec_bin(kn,bn,nx,ny,kmap,bp,dbin_max)

def binfilter(nx,ny,D,bp,bn=1):
  """
  The multipole bin binary-mask given on the 2D Fourier grids so that M = 1 inside the multipole bin and 0 otherwise

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*2*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :bp[*bn+1*] (*double*): Multipole bin edges

  Returns:
    :bf[*bn,nx,ny*] (*double*): The multipole bin binary-mask for each multipol bin
  Usage:
    :bf = flatsky.bispec.binfilter(nx,ny,D,bp,bn):
  """
  bn = len(bp) - 1
  return libflatsky.bispec.binfilter(nx,ny,D,bp,bn)

def bispec_norm_1d(nx,ny,D,bfs,bn=1):
  """
  Normalization of the 1D binned bispectrum estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*2*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :bfs[*3,bn,nx,ny*] (*double*): Multipole bin mask on 2D grids obtained from the binfilter function

  Returns:
    :bnorm[*bn*] (*double*): Normalization of 1D binned bispectrum at each multipole bin
  Usage:
    :bnorm = flatsky.bispec.bispec_norm_1d(nx,ny,D,bfs,bn):
  """
  bn = len(bfs[0,:,0,0])
  return libflatsky.bispec.bispec_norm_1d(nx,ny,D,bfs,bn)

def bispec_bin_1d(nx,ny,D,bfs,bnorm,alm,bn=1):
  """
  1D binned bispectrum estimator

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*2*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :bfs[*3,bn,nx,ny*] (*double*): Multipole bin mask on 2D grids obtained from the binfilter function
    :bnorm[*bn*] (*double*): Normalization of 1D binned bispectrum at each multipole bin
    :alm[*3,nx,ny*] (*dcmplx*): Fourier modes for each leg

  Returns:
    :bispec[*bn*] (*double*): 1D binned bispectrum at each multipole bin
  Usage:
    :bispec = flatsky.bispec.bispec_bin_1d(nx,ny,D,bfs,bnorm,alm,bn):
  """
  bn = len(bfs[0,:,0,0])
  return libflatsky.bispec.bispec_bin_1d(nx,ny,D,bfs,bnorm,alm,bn)

