import libflatsky

def dft1d(map0,nx,ny,D,trans):
  """
  Usage:
    - e.g., npix,map1 = flatsky.ffttools.dft1d(map0,nx,ny,D,trans)
  """
  return flatsky.ffttools.dft1d(map0,nx,ny,D,trans)

def dft2d(map0,nx,ny,D,trans):
  """
  Usage:
    - e.g., map1 = flatsky.ffttools.dft2d(map0,nx,ny,D,trans)
  """
  return flatsky.ffttools.dft2d(map0,nx,ny,D,trans)

def dft2dr(map0,nx,ny,D,trans):
  """
  Usage:
    - e.g., map1 = flatsky.ffttools.dft2dr(map0,nx,ny,D,trans)
  """
  return flatsky.ffttools.dft2dr(map0,nx,ny,D,trans)

def dft2drc(map0,nx,ny,D,trans):
  """
  Usage:
    - e.g., map1 = flatsky.ffttools.dft2drc(map0,nx,ny,D,trans)
  """
  return flatsky.ffttools.dft2drc(map0,nx,ny,D,trans)

def dft2dcr(map0,nx,ny,D,trans):
  """
  Usage:
    - e.g., map1 = flatsky.ffttools.dft2dcr(map0,nx,ny,D,trans)
  """
  return flatsky.ffttools.dft2dcr(map0,nx,ny,D,trans)

def eb_separate(nx,ny,D,QU,W=0,Wd=0):
  """
   Transform to Fourier Space
   add corrections
  Usage:
    - e.g., EB = flatsky.ffttools.eb_separate(nx,ny,D,QU,W,Wd)
  """
  if W==0: W= 0
  if Wd==0: Wd= 0
  return flatsky.ffttools.eb_separate(nx,ny,D,QU,W,Wd)

