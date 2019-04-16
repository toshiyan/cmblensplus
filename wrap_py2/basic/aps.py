import libbasic

def binning(bn,eL,spc=0):
  """
  Usage:
    - e.g., bp,bc = basic.aps.binning(bn,eL,spc)
  """
  if spc==0: spc= ''
  return basic.aps.binning(bn,eL,spc)

def read_cambcls(f,lmin,lmax,numcls,bb=0,raw=0):
  """
  Usage:
    - e.g., cl = basic.aps.read_cambcls(f,lmin,lmax,numcls,bb,raw)
  """
  if bb==0: bb= 0
  if raw==0: raw= 0
  return basic.aps.read_cambcls(f,lmin,lmax,numcls,bb,raw)

