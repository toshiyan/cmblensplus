import libanalytic

def alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe=0,gln=0,gle=0,lxcut=0):
  """
  Usage:
    - e.g., Ag,Ac = analytic.flat.alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe,gln,gle,lxcut)
  """
  if qe==0: qe= 'lensing'
  if gln==0: gln= 100
  if gle==0: gle= 1e-14
  if lxcut==0: lxcut= 0
  return analytic.flat.alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe,gln,gle,lxcut)

def bbxy(lmax,rlmin,rlmax,XX,YY,weight=0,gln=0,gle=0):
  """
  Usage:
    - e.g., BB = analytic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)
  """
  if weight==0: weight= 'lensing'
  if gln==0: gln= 100
  if gle==0: gle= 1e-14
  return analytic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)

