import libflatsky

def qtt(nx,ny,D,rL,T1,T2):
  """
  Usage:
    :slm = flatsky.rec_src.qtt(nx,ny,D,rL,T1,T2):
  """
  return libflatsky.rec_src.qtt(nx,ny,D,rL,T1,T2)

def qte(nx,ny,D,rL,T,E):
  """
  Usage:
    :slm = flatsky.rec_src.qte(nx,ny,D,rL,T,E):
  """
  return libflatsky.rec_src.qte(nx,ny,D,rL,T,E)

def qtb(nx,ny,D,rL,T,B):
  """
  Usage:
    :slm = flatsky.rec_src.qtb(nx,ny,D,rL,T,B):
  """
  return libflatsky.rec_src.qtb(nx,ny,D,rL,T,B)

def qee(nx,ny,D,rL,E1,E2):
  """
  Usage:
    :slm = flatsky.rec_src.qee(nx,ny,D,rL,E1,E2):
  """
  return libflatsky.rec_src.qee(nx,ny,D,rL,E1,E2)

def qeb(nx,ny,D,rL,E,B):
  """
  Usage:
    :slm = flatsky.rec_src.qeb(nx,ny,D,rL,E,B):
  """
  return libflatsky.rec_src.qeb(nx,ny,D,rL,E,B)

