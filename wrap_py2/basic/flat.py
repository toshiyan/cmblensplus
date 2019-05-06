import libbasic

def alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe=None,gln=None,gle=None,lxcut=None):
  """
  Usage:
    - e.g., Ag,Ac = basic.flat.alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe,gln,gle,lxcut)
  """
  if qe is None: qe= 'lensing'
  if gln is None: gln= 100
  if gle is None: gle= 1e-14
  if lxcut is None: lxcut= 0
  return libbasic.flat.alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe,gln,gle,lxcut)

def bbxy(lmax,rlmin,rlmax,XX,YY,weight=None,gln=None,gle=None):
  """
  Usage:
    - e.g., BB = basic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)
  """
  if weight is None: weight= 'lensing'
  if gln is None: gln= 100
  if gle is None: gle= 1e-14
  return libbasic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)

