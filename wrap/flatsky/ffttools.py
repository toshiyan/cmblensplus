import libflatsky

def dft1d(map0,nx,ny,npix,D,trans,map1):
  """
  DFT for 1D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :npix (*int*): Total number of grids (npix=nx*ny)
    :trans (*int*): 1 (map to Fourier) or -1 (Fourier to map)
    :D[*2*] (*double*): Side length (x and y) of map
    :map0[*pix*] (*dcmplx*): Data on 2D grid to be transformed, with bounds (npix)

  Returns:
    :map1[*pix*] (*dcmplx*): Transformed data on 2D grid, with bounds (npix)

  Usage:
    :map1(npix) = flatsky.ffttools.dft1d(map0,nx,ny,npix,D,trans,map1):
  """
  return libflatsky.ffttools.dft1d(map0,nx,ny,npix,D,trans,map1)

def dft2d(map0,nx,ny,D,trans):
  """
  DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :trans (*int*): 1 (map to Fourier) or -1 (Fourier to map)
    :D[*2*] (*double*): Side length (x and y) of map
    :map0[*x,y*] (*dcmplx*): Data on 2D grid to be transformed, with bounds (nx,ny)

  Returns:
    :map1[*x,y*] (*dcmplx*): Transformed data on 2D grid, with bounds (nx,ny)

  Usage:
    :map1 = flatsky.ffttools.dft2d(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2d(map0,nx,ny,D,trans)

def dft2dr(map0,nx,ny,D,trans):
  """
  DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :trans (*int*): 1 (map to Fourier) or -1 (Fourier to map)
    :D[*2*] (*double*): Side length (x and y) of map
    :map0[*x,y*] (*double*): Data on 2D grid to be transformed, with bounds (nx,ny)

  Returns:
    :map1[*x,y*] (*double*): Transformed data on 2D grid, with bounds (nx,ny)

  Usage:
    :map1 = flatsky.ffttools.dft2dr(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2dr(map0,nx,ny,D,trans)

def dft2drc(map0,nx,ny,D,trans):
  """
  DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :trans (*int*): 1 (map to Fourier) or -1 (Fourier to map)
    :D[*2*] (*double*): Side length (x and y) of map
    :map0[*x,y*] (*double*): Data on 2D grid to be transformed, with bounds (nx,ny)

  Returns:
    :map1[*x,y*] (*dcmplx*): Transformed data on 2D grid, with bounds (nx,ny)

  Usage:
    :map1 = flatsky.ffttools.dft2drc(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2drc(map0,nx,ny,D,trans)

def dft2dcr(map0,nx,ny,D,trans):
  """
  DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :trans (*int*): 1 (map to Fourier) or -1 (Fourier to map)
    :D[*2*] (*double*): Side length (x and y) of map
    :map0[*x,y*] (*dcmplx*): Data on 2D grid to be transformed, with bounds (nx,ny)

  Returns:
    :map1[*x,y*] (*double*): Transformed data on 2D grid, with bounds (nx,ny)

  Usage:
    :map1 = flatsky.ffttools.dft2dcr(map0,nx,ny,D,trans):
  """
  return libflatsky.ffttools.dft2dcr(map0,nx,ny,D,trans)

def dft2dpol(nx,ny,D,Q,U):
  """
  Spin-2 DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :D[*2*] (*double*): Side length (x and y) of map
    :Q[*x,y*] (doub;e): Q on 2D grid to be transformed, with bounds (nx,ny)
    :U[*x,y*] (doub;e): U on 2D grid to be transformed, with bounds (nx,ny)

  Returns:
    :E[*x,y*] (*dcmplx*): E on 2D Fourier grid, with bounds (nx,ny)
    :B[*x,y*] (*dcmplx*): B on 2D Fourier grid, with bounds (nx,ny)

  Usage:
    :E,B = flatsky.ffttools.dft2dpol(nx,ny,D,Q,U):
  """
  return libflatsky.ffttools.dft2dpol(nx,ny,D,Q,U)

def idft2dpol(nx,ny,D,E,B):
  """
  Spin-2 Inverse DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :D[*2*] (*double*): Side length (x and y) of map
    :E[*x,y*] (*dcmplx*): E on 2D Fourier grid, with bounds (nx,ny)
    :B[*x,y*] (*dcmplx*): B on 2D Fourier grid, with bounds (nx,ny)

  Returns:
    :Q[*x,y*] (doub;e): Q on 2D grid, with bounds (nx,ny)
    :U[*x,y*] (doub;e): U on 2D grid, with bounds (nx,ny)

  Usage:
    :Q,U = flatsky.ffttools.idft2dpol(nx,ny,D,E,B):
  """
  return libflatsky.ffttools.idft2dpol(nx,ny,D,E,B)

def eb_separate(nx,ny,D,QU,W,Wd=None):
  """
  Compute Smith's pure EB estimator in flatsky

  Args:
    :nx, ny (*int*): Number of x/Lx and y/Ly grids
    :D[*2*] (*double*): Side length (x and y) of map
    :W[*x,y*] (*double*): Window function, with bounds (nx,ny)
    :QU[*x,y,2*] (*double*): unmasked Q and U maps, with bounds (nx,ny,2)

  Args(optional):
    :Wd[*5,x,y*] (*double*): Precomputed window function derivaives, dW/dx, dW/dw, d^2W/dx^2, d^2W/dxdy, d^2W/dy^2, with bounds (5,nx,ny)

  Returns:
    :EB[*2,x,y*] (*dcmplx*): E and B modes in 2D Fourier grid, with bounds (2,nx,ny)

  Usage:
    :EB = flatsky.ffttools.eb_separate(nx,ny,D,QU,W,Wd):
  """
  if Wd is None: Wd= np.zeros((5,nx,ny))
  return libflatsky.ffttools.eb_separate(nx,ny,D,QU,W,Wd)

