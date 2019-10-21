import libbasic

def alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,gln= 100,gle= 1e-14,lxcut= 0):
  """
  Usage:
    :Ag,Ac = basic.flat.alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,gln,gle,lxcut):
  """
  return libbasic.flat.alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,gln,gle,lxcut)

def alxy_assym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,gln= 100,gle= 1e-14,lxcut= 0):
  """
  Usage:
    :Ag,Ac = basic.flat.alxy_assym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,gln,gle,lxcut):
  """
  return libbasic.flat.alxy_assym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,gln,gle,lxcut)

def bbxy(lmax,rlmin,rlmax,XX,YY,weight= 'lensing',gln= 100,gle= 1e-14):
  """
  Usage:
    :BB = basic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle):
  """
  return libbasic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)

