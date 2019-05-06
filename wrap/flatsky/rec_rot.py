import libflatsky

def qte(nx,ny,D,rL,fC,T,E):
  """
  Reconstructing anisotropic pol. rot. angles from TE quadratic estimator
  Usage:
    :alm = flatsky.rec_rot.qte(nx,ny,D,rL,fC,T,E):
  """
  return libflatsky.rec_rot.qte(nx,ny,D,rL,fC,T,E)

def qtb(nx,ny,D,rL,fC,T,B):
  """
  Reconstructing anisotropic pol. rot. angles from TB quadratic estimator
  Usage:
    :alm = flatsky.rec_rot.qtb(nx,ny,D,rL,fC,T,B):
  """
  return libflatsky.rec_rot.qtb(nx,ny,D,rL,fC,T,B)

def qee(nx,ny,D,rL,fC,E1,E2):
  """
  Reconstructing anisotropic pol. rot. angles from EE quadratic estimator
  Usage:
    :alm = flatsky.rec_rot.qee(nx,ny,D,rL,fC,E1,E2):
  """
  return libflatsky.rec_rot.qee(nx,ny,D,rL,fC,E1,E2)

def qeb(nx,ny,D,rL,EE,E,B=None,BB=None):
  """
  Reconstructing anisotropic pol. rot. angles from EB quadratic estimator
  Usage:
    :alm = flatsky.rec_rot.qeb(nx,ny,D,rL,EE,E,B,BB):
  """
  if B is None: B= ''
  if BB is None: BB= ''
  return libflatsky.rec_rot.qeb(nx,ny,D,rL,EE,E,B,BB)

