import libbasic

def binning(bn,eL,spc=''):
  """
 Return multipole-bin edges and centers
 
 Args:
    :bn (*int*): number of bins
    :eL[*2*] (*int*): bin edges

 Args(optional):
    :spc (*str*): bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing

 Returns:
    :bp (*double*): bin edges, with bounds (0:bn)
    :bc (*double*): bin centers, with bounds (bn)

  Usage:
    :bp,bc = basic.aps.binning(bn,eL,spc):
  """
  return libbasic.aps.binning(bn,eL,spc)

def read_cambcls(f,lmin,lmax,numcls,bb=False,lsq=True):
  """
  Return CMB cls from CAMB output files

  Args:
    :f (*str*): Filename
    :lmin (*int*): Minimum multipole of the output Cl
    :lmax (*int*): Maximum multipole of the output Cl
    :numcls (*int*): Number of cls to be read

  Args(Optional):
    :bb (*bool*): Filetype. The data should be TT, EE, TE, dd, Td (,Ed) if bb = False (default), and TT, EE, BB, TE if bb = True.
    :lsq (*bool*): The cls is divided by l(l+1)/2pi if lsq = True (default)

  Returns:
    :cl[*numcls,l*] (*double*): Angular power spectra, with bountd (numcls,0:lmax)

  Usage:
    :cl = basic.aps.read_cambcls(f,lmin,lmax,numcls,bb,lsq):
  """
  return libbasic.aps.read_cambcls(f,lmin,lmax,numcls,bb,lsq)

def map_vars(lmax,cl):
  """
  Variance of a map and its derivative
 
  Args:
    :lmax (*int*): Maximum multipole of the input cl
    :cl[*l*] (*double*): Angular power spectrum, with bounds (0:lmax)

  Returns:
    :sigma[*2*] (*double*): Variance of the map and of map derivative, with bounds (2)
  
  Usage:
    :sigma = basic.aps.map_vars(lmax,cl):
  """
  return libbasic.aps.map_vars(lmax,cl)

def cl2bcl(bn,lmax,cl,lmin=1,spc=''):
  """
  From unbinned to binned angular power spectrum

  Args:
    :bn (*int*): number of multipole bins
    :lmax (*int*): maximum multipole of the input angular power spectrum
    :cl[*l*] (*double*): angular power spectrum, with bounds (0:lmax)

 Args(optional):
    :lmin (*int*): minimum multipole
    :spc (*str*): bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing

  Returns:
    :cb[*bin*] (*double*): auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)

  Usage:
    :cb = basic.aps.cl2bcl(bn,lmax,cl,lmin,spc):
  """
  return libbasic.aps.cl2bcl(bn,lmax,cl,lmin,spc)

