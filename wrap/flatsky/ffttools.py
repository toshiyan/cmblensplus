import libflatsky

def dft1d(map0,nx,ny,D,trans):
  """
  Usage:
    :npix,map1 = flatsky.ffttools.dft1d(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft1d(map0,nx,ny,D,trans)

def dft2d(map0,nx,ny,D,trans):
  """
  Usage:
    :map1 = flatsky.ffttools.dft2d(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2d(map0,nx,ny,D,trans)

def dft2dr(map0,nx,ny,D,trans):
  """
  Usage:
    :map1 = flatsky.ffttools.dft2dr(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2dr(map0,nx,ny,D,trans)

def dft2drc(map0,nx,ny,D,trans):
  """
  Usage:
    :map1 = flatsky.ffttools.dft2drc(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2drc(map0,nx,ny,D,trans)

def dft2dcr(map0,nx,ny,D,trans):
  """
  Usage:
    :map1 = flatsky.ffttools.dft2dcr(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2dcr(map0,nx,ny,D,trans)

def eb_separate(nx,ny,D,QU,W=None,Wd=None):
  """
   Transform to Fourier Space
   add corrections
  Usage:
    :EB = flatsky.ffttools.eb_separate(nx,ny,D,QU,W,Wd):
  """
  if W is None: W= 0
  if Wd is None: Wd= 0
  return libflatsky.ffttools.eb_separate(nx,ny,D,QU,W,Wd)

