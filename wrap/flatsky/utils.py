import libflatsky
import numpy

def map2alm(nx,ny,D,map):
  """
  DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of x and y grids
    :D[*2*] (*double*): Side length (x and y) of map
    :map[*x,y*] (*double*): Map on 2D grid with bounds (nx,ny)

  Returns:
    :alm[*x,y*] (*dcmplx*): Fourier modes on 2D grid, with bounds (nx,ny)

  Usage:
    :alm = flatsky.utils.map2alm(nx,ny,D,map):
  """
  return libflatsky.utils.map2alm(nx,ny,D,map)

def alm2map(nx,ny,D,alm):
  """
  DFT for 2D array. 

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*2*] (*double*): Side length (x and y) of map
    :alm[*x,y*] (*dcmplx*): Fourier modes on 2D grid to be transformed, with bounds (nx,ny)

  Returns:
    :map[*x,y*] (*double*): Map on 2D grid, with bounds (nx,ny)

  Usage:
    :map = flatsky.utils.alm2map(nx,ny,D,alm):
  """
  return libflatsky.utils.alm2map(nx,ny,D,alm)

def el2d(nx,ny,D):
  """
  Return absolute value of multipole in 2D grids
 
  Args:
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
  
  Returns:
    :els[*nx,ny*] (*double*): absolute value of Fourier mode, (Lx**2+Ly**2)**0.5, with bounds (nx,ny)

  Usage:
    :els = flatsky.utils.el2d(nx,ny,D):
  """
  return libflatsky.utils.el2d(nx,ny,D)

def elarrays(nx,ny,D):
  """
  Return Lx, Ly, absolute value of multipole, and its inverse in 2D grids
 
  Args:
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
  
  Returns:
    :elx[*nx,ny*] (*double*): Lx, with bounds (nx,ny)
    :ely[*nx,ny*] (*double*): Ly, with bounds (nx,ny)
    :els[*nx,ny*] (*double*): absolute value of Fourier mode, (Lx**2+Ly**2)**0.5, with bounds (nx,ny)
    :eli[*nx,ny*] (*double*): inverse of els, with bounds (nx,ny)

  Usage:
    :elx,ely,els,eli = flatsky.utils.elarrays(nx,ny,D):
  """
  return libflatsky.utils.elarrays(nx,ny,D)

def elmask(nx,ny,D,lmin= 0,lmax= 1000,lxcut= 0,lycut= 0):
  """
  Return mask in 2D Fourier space. The lmask is unity at lmin<=|L|<=lmax, |Lx|>=lxcut, |Ly|>=lycut, and otherwize zero. 

  Args: 
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)

  Args(optional):
    :lmin/lmax (*int*): Minimum/Maximum of multipoles
    :lxcut/lycut (*int*): Remove |Lx|<lxcut / |Ly|<lycut cut of multipoles
  
  Returns:
    :lmask[*nx,ny*] (*double*): Mask, with bounds (nx,ny)

  Usage:
    :lmask = flatsky.utils.elmask(nx,ny,D,lmin,lmax,lxcut,lycut):
  """
  return libflatsky.utils.elmask(nx,ny,D,lmin,lmax,lxcut,lycut)

def alm2bcl(bn,oL,nx,ny,D,alm1,spc='',alm2=None):
  """
  Compute angular power spectrum from Fourier modes, with multipole binning
 
  Args:
    :bn (*int*): number of multipole bin
    :oL[*2*] (*int*): minimum and maximum multipoles of the output cl
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :alm1[*nx,ny*] (*dcmplx*): Fourier mode, with bounds (nx,ny)
 
  Args(optional):
    :alm2[*nx,ny*] (*dcmplx*): Fourier mode, with bounds (nx,ny), default to None
    :spc (*str*): type of multipole binning, i.e., linear spacing (spc='', default), or log spacing (spc='log')

  Returns:
    :Cb[*bin*] (*double*): angular power spectrum with multipole binning, with bounds (bn)

  Usage:
    :Cb = flatsky.utils.alm2bcl(bn,oL,nx,ny,D,alm1,alm2,spc):
  """
  if alm2 is None: alm2= alm1
  return libflatsky.utils.alm2bcl(bn,oL,nx,ny,D,alm1,alm2,spc)

def c2d2bcl(nx,ny,D,c2d,bn,oL,spc=''):
  """
  Return 1D angular power spectrum with multipole binning from a 2D power spectrum

  Args:
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :c2d[*nx,ny*] (*double*): 2D power spectrum, with bounds (nx,ny)
    :bn (*int*): number of multipole bin
    :oL[*2*] (*int*): minimum and maximum multipoles of the output cl
    
  Args(optional):
    :spc (*str*): type of multipole binning, i.e., linear spacing (spc='', default), or log spacing (spc='log')

  Returns:
    :Cb[*bin*] (*double*): angular power spectrum with multipole binning, with bounds (bn)

  Usage:
    :Cb = flatsky.utils.c2d2bcl(nx,ny,D,c2d,bn,oL,spc):
  """
  return libflatsky.utils.c2d2bcl(nx,ny,D,c2d,bn,oL,spc)

def cl2c2d(nx,ny,D,lmin,lmax,Cl,method='linear'):
  """
  Assign values of 1D angular power spectrum on to 2D grid with linear interpolation

  Args: 
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :lmin (*int*): minimum multipole of cl to be interpolated
    :lmax (*int*): maximum multipole of cl to be interpolated
    :Cl[*l*] (*double*): 1D power spectrum, with bounds (0:lmax)

  Args(optional):
    :method (*str*): type of interpolation method, i.e., linear interpolation (method='linear', default), or step (method='step')

  Returns:
    :c2d[*nx,ny*] (*double*): 2D power spectrum, with bounds (nx,ny)
 
  Usage:
    :c2d = flatsky.utils.cl2c2d(nx,ny,D,lmin,lmax,Cl,method):
  """
  return libflatsky.utils.cl2c2d(nx,ny,D,lmin,lmax,Cl,method)

def cb2c2d(bn,bc,nx,ny,D,lmin,lmax,Cb,method=''):
  """
  Assign values of 1D angular power spectrum on to 2D grid with linear interpolation

  Args: 
    :bn (*int*): number of multipole bins
    :bc[*bin*] (*double*): multipole bin center, with bounds (bn)
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :lmin (*int*): minimum multipole of cl to be interpolated
    :lmax (*int*): maximum multipole of cl to be interpolated
    :Cb[*bin*] (*double*): 1D power spectrum with multipole binning, with bounds (bn)

  Args(optional):
    :method (*str*): interpolation method from binned to unbinned angular spectrum, i.e., spline ('', default), or linear ('linear') interpolation

  Returns:
    :c2d[*nx,ny*] (*double*): 2D power spectrum, with bounds (nx,ny)
 
  Usage:
    :C2d = flatsky.utils.cb2c2d(bn,bc,nx,ny,D,lmin,lmax,Cb,method):
  """
  return libflatsky.utils.cb2c2d(bn,bc,nx,ny,D,lmin,lmax,Cb,method)

def gauss1alm(nx,ny,D,lmin,lmax,Cl):
  """
  Generate random gaussian fields in 2D Fourier space for a given isotropic spectrum, satisfying a^*_l = a_{-l}

  Args:
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :lmin (*int*): minimum multipole of cl to be interpolated
    :lmax (*int*): maximum multipole of cl to be interpolated
    :Cl[*l*] (*double*): 1D power spectrum, with bounds (0:lmax)

  Returns:
    :alm[*lx,ly*] (*dcmplx*): random gaussian fields on 2D Fourier plane, with bounds (nx,ny)

  Usage:
    :alm = flatsky.utils.gauss1alm(nx,ny,D,lmin,lmax,Cl):
  """
  return libflatsky.utils.gauss1alm(nx,ny,D,lmin,lmax,Cl)

def gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE):
  """
  Generate two correlated random gaussian fields in 2D Fourier space for a given isotropic spectrum

  Args:
    :nx, ny (*int*): number of Lx and Ly grids
    :D[*xy*] (*double*): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
    :lmin (*int*): minimum multipole of cl to be interpolated
    :lmax (*int*): maximum multipole of cl to be interpolated
    :TT[*l*] (*double*): the 1st 1D power spectrum, with bounds (0:lmax)
    :TE[*l*] (*double*): the cross 1D power spectrum, with bounds (0:lmax)
    :EE[*l*] (*double*): the 2nd 1D power spectrum, with bounds (0:lmax)

  Returns:
    :tlm[*lx,ly*] (*dcmplx*): the 1st random gaussian fields on 2D Fourier plane, with bounds (nx,ny)
    :elm[*lx,ly*] (*dcmplx*): the 2nd random gaussian fields on 2D Fourier plane, with bounds (nx,ny)

  Usage:
    :tlm,elm = flatsky.utils.gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE):
  """
  return libflatsky.utils.gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE)

def window_sin(nx,ny,D,ap= 1,cut= 1):
  """
  Return a sin window function.

  Args:
    :nx, ny (*int*): Number of Lx and Ly grids
    :D[*xy*] (*double*): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)

  Args(Optional):
    :ap (*double*): Apodization parameter defined by apodized-range = (1-ap) x (cut)mapsize, from 0 (full apodization) to 1 (no apodization). Default to 1.
    :cut (*double*): Map cut scale defined by cutmapsize = cut x mapsize, from 0 (full cut) to 1 (no cut). Default to 1.

  Return:
    :W[*x,y*] (*double*): Window function, with bounds (nx,ny)

  Usage:
    :W = flatsky.utils.window_sin(nx,ny,D,ap,cut):
  """
  return libflatsky.utils.window_sin(nx,ny,D,ap,cut)

def window_norm(nx,ny,wind,num):
  """
  Usage:
    :wn = flatsky.utils.window_norm(nx,ny,wind,num):
  """
  return libflatsky.utils.window_norm(nx,ny,wind,num)

def window_norm_x(nx,ny,W1,W2,num):
  """
  Usage:
    :Wn = flatsky.utils.window_norm_x(nx,ny,W1,W2,num):
  """
  return libflatsky.utils.window_norm_x(nx,ny,W1,W2,num)

def rotation(nx,ny,rot,QU,rtype):
  """
  Usage:
    :rQU = flatsky.utils.rotation(nx,ny,rot,QU,rtype):
  """
  return libflatsky.utils.rotation(nx,ny,rot,QU,rtype)

def get_angle(nx,ny,D):
  """
  Usage:
    :theta,phi = flatsky.utils.get_angle(nx,ny,D):
  """
  return libflatsky.utils.get_angle(nx,ny,D)

def cutmap(ox,oy,cx,cy,omap):
  """
  Usage:
    :cmap = flatsky.utils.cutmap(ox,oy,cx,cy,omap):
  """
  return libflatsky.utils.cutmap(ox,oy,cx,cy,omap)

