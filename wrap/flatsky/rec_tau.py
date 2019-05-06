import libflatsky

def qtt(nx,ny,D,rL,fC,T1,T2):
  """
  Usage:
    :tlm = flatsky.rec_tau.qtt(nx,ny,D,rL,fC,T1,T2):
  """
  return libflatsky.rec_tau.qtt(nx,ny,D,rL,fC,T1,T2)

def qte(nx,ny,D,rL,fC,T,E):
  """
  Usage:
    :tlm = flatsky.rec_tau.qte(nx,ny,D,rL,fC,T,E):
  """
  return libflatsky.rec_tau.qte(nx,ny,D,rL,fC,T,E)

def qtb(nx,ny,D,rL,fC,T,B):
  """
  Usage:
    :tlm = flatsky.rec_tau.qtb(nx,ny,D,rL,fC,T,B):
  """
  return libflatsky.rec_tau.qtb(nx,ny,D,rL,fC,T,B)

def qee(nx,ny,D,rL,fC,E1,E2):
  """
  Usage:
    :tlm = flatsky.rec_tau.qee(nx,ny,D,rL,fC,E1,E2):
  """
  return libflatsky.rec_tau.qee(nx,ny,D,rL,fC,E1,E2)

def qeb(nx,ny,D,rL,fC,E,B):
  """
  Usage:
    :tlm = flatsky.rec_tau.qeb(nx,ny,D,rL,fC,E,B):
  """
  return libflatsky.rec_tau.qeb(nx,ny,D,rL,fC,E,B)

