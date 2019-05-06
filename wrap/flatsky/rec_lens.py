import libflatsky

def qtt(nx,ny,D,rL,fC,T1,T2,gtype=None):
  """
  Usage:
    :glm,clm = flatsky.rec_lens.qtt(nx,ny,D,rL,fC,T1,T2,gtype):
  """
  if gtype is None: gtype= ''
  return libflatsky.rec_lens.qtt(nx,ny,D,rL,fC,T1,T2,gtype)

def qte(nx,ny,D,rL,fC,T,E,gtype=None):
  """
  Usage:
    :glm,clm = flatsky.rec_lens.qte(nx,ny,D,rL,fC,T,E,gtype):
  """
  if gtype is None: gtype= ''
  return libflatsky.rec_lens.qte(nx,ny,D,rL,fC,T,E,gtype)

def qtb(nx,ny,D,rL,fC,T,B,gtype=None):
  """
  Usage:
    :glm,clm = flatsky.rec_lens.qtb(nx,ny,D,rL,fC,T,B,gtype):
  """
  if gtype is None: gtype= ''
  return libflatsky.rec_lens.qtb(nx,ny,D,rL,fC,T,B,gtype)

def qee(nx,ny,D,rL,fC,E1,E2,gtype=None):
  """
  Usage:
    :glm,clm = flatsky.rec_lens.qee(nx,ny,D,rL,fC,E1,E2,gtype):
  """
  if gtype is None: gtype= ''
  return libflatsky.rec_lens.qee(nx,ny,D,rL,fC,E1,E2,gtype)

def qeb(nx,ny,D,rL,fC,E,B,gtype=None):
  """
  Usage:
    :glm,clm = flatsky.rec_lens.qeb(nx,ny,D,rL,fC,E,B,gtype):
  """
  if gtype is None: gtype= ''
  return libflatsky.rec_lens.qeb(nx,ny,D,rL,fC,E,B,gtype)

def qbb(nx,ny,D,rL,fC,B1,B2,gtype=None):
  """
  Usage:
    :glm,clm = flatsky.rec_lens.qbb(nx,ny,D,rL,fC,B1,B2,gtype):
  """
  if gtype is None: gtype= ''
  return libflatsky.rec_lens.qbb(nx,ny,D,rL,fC,B1,B2,gtype)

