import libflatsky

def qeb(nx,ny,D,rL,IE,IB,EE,eL,BB=None):
  """
  Return normalization of EB quadratic estimator for anisotropic pol. rot. angles
  Usage:
    :Aa = flatsky.norm_rot.qeb(nx,ny,D,rL,IE,IB,EE,eL,BB):
  """
  if BB is None: BB= ''
  return libflatsky.norm_rot.qeb(nx,ny,D,rL,IE,IB,EE,eL,BB)

