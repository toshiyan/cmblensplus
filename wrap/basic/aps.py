import libbasic

def binning(bn,eL,spc=None):
  """
 Return multipole-bin edges and centers
 
 Args:
    - bn (*int*): number of bins
    - eL[*2*] (*int*): bin edges

 Args(optional):
    - spc (*str*): bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing

 Returns:
    - bp (*double*): bin edges, with bounds (0:bn)
    - bc (*double*): bin centers, with bounds (bn)

  Usage:
    :bp,bc = basic.aps.binning(bn,eL,spc):
  """
  if spc is None: spc= ''
  return libbasic.aps.binning(bn,eL,spc)

def read_cambcls(f,lmin,lmax,numcls,bb=None,raw=None):
  """
  Usage:
    :cl = basic.aps.read_cambcls(f,lmin,lmax,numcls,bb,raw):
  """
  if bb is None: bb= 0
  if raw is None: raw= 0
  return libbasic.aps.read_cambcls(f,lmin,lmax,numcls,bb,raw)

def map_vars(lmax,cl):
  """
  Usage:
    :sigma = basic.aps.map_vars(lmax,cl):
  """
  return libbasic.aps.map_vars(lmax,cl)

def cl2bcl(bn,lmax,cl,spc=None):
  """
  From unbinned to binned angular power spectrum

  Args:
    - bn (*int*): number of multipole bins
    - lmax (*int*): maximum multipole of the input angular power spectrum
    - cl[*l*] (*double*): angular power spectrum, with bounds (0:lmax)

 Args(optional):
    - spc (*str*): bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing

  Returns:
    - cb[*bin*] (*double*): auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)

  Usage:
    :cb = basic.aps.cl2bcl(bn,lmax,cl,spc):
  """
  if spc is None: spc= ''
  return libbasic.aps.cl2bcl(bn,lmax,cl,spc)

