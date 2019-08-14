import libbasic

def alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe= 'lensing',gln= 100,gle= 1e-14,lxcut= 0):
  """
  Usage:
    :Ag,Ac = basic.flat.alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe,gln,gle,lxcut):
  """
  return libbasic.flat.alxy(est,lmax,rlmin,rlmax,fC,W1,W2,qe,gln,gle,lxcut)

def bbxy(lmax,rlmin,rlmax,XX,YY,weight= 'lensing',gln= 100,gle= 1e-14):
  """
  Usage:
    :BB = basic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle):
  """
  return libbasic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)

