import libflatsky

def el2d(nx,ny,D):
  """
  Usage:
    - e.g., els = flatsky.utils.el2d(nx,ny,D)
  """
  return flatsky.utils.el2d(nx,ny,D)

def elarrays(nx,ny,D):
  """
  Usage:
    - e.g., elx,ely,els,eli = flatsky.utils.elarrays(nx,ny,D)
  """
  return flatsky.utils.elarrays(nx,ny,D)

def alm2bcl(bn,oL,nx,ny,D,alm1,alm2=0,spc=0):
  """
   2D power spectrum
   to 1D power spectrum
  Usage:
    - e.g., Cb = flatsky.utils.alm2bcl(bn,oL,nx,ny,D,alm1,alm2,spc)
  """
  if alm2==0: alm2= alm1
  if spc==0: spc=''
  return flatsky.utils.alm2bcl(bn,oL,nx,ny,D,alm1,alm2,spc)

def c2d2bcl(nx,ny,D,Cl,bn,oL,spc=0):
  """
  Usage:
    - e.g., Cb = flatsky.utils.c2d2bcl(nx,ny,D,Cl,bn,oL,spc)
  """
  if spc==0: spc=''
  return flatsky.utils.c2d2bcl(nx,ny,D,Cl,bn,oL,spc)

def cl2c2d(nx,ny,D,lmin,lmax,Cl):
  """
 Transform Cl to Cl2D with linear interpolation
  Usage:
    - e.g., c2d = flatsky.utils.cl2c2d(nx,ny,D,lmin,lmax,Cl)
  """
  return flatsky.utils.cl2c2d(nx,ny,D,lmin,lmax,Cl)

def cb2c2d(bn,bc,nx,ny,D,eL,Cb,method0=0,bp=0):
  """
  Usage:
    - e.g., C2d = flatsky.utils.cb2c2d(bn,bc,nx,ny,D,eL,Cb,method0,bp)
  """
  if method0==0: method0=''
  if bp==0: bp=''
  return flatsky.utils.cb2c2d(bn,bc,nx,ny,D,eL,Cb,method0,bp)

def gauss1alm(nx,ny,D,lmin,lmax,Cl):
  """
   make cl on 2d grid
   alm=0 if i=1 or j=1 for symmetry
   center: (ic,jc) = (nx/2+1, ny/2+1)
   maximum nn
   check
  Usage:
    - e.g., alm = flatsky.utils.gauss1alm(nx,ny,D,lmin,lmax,Cl)
  """
  return flatsky.utils.gauss1alm(nx,ny,D,lmin,lmax,Cl)

def gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE):
  """
  Usage:
    - e.g., tlm,elm = flatsky.utils.gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE)
  """
  return flatsky.utils.gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE)

def window_sin(nx,ny,D,ap=0,cut=0):
  """
  Usage:
    - e.g., W = flatsky.utils.window_sin(nx,ny,D,ap,cut)
  """
  if ap==0: ap= 1
  if cut==0: cut= 1
  return flatsky.utils.window_sin(nx,ny,D,ap,cut)

def window_norm(nx,ny,numn):
  """
  Usage:
    - e.g., W,Wn = flatsky.utils.window_norm(nx,ny,numn)
  """
  return flatsky.utils.window_norm(nx,ny,numn)

def window_norm_x(nx,ny,W1,W2,num):
  """
  Usage:
    - e.g., Wn = flatsky.utils.window_norm_x(nx,ny,W1,W2,num)
  """
  return flatsky.utils.window_norm_x(nx,ny,W1,W2,num)

def rotation(nx,ny,rot,rtype):
  """
  Usage:
    - e.g., QU,rQU = flatsky.utils.rotation(nx,ny,rot,rtype)
  """
  return flatsky.utils.rotation(nx,ny,rot,rtype)

def get_angle(nx,ny,D):
  """
  Usage:
    - e.g., theta,phi = flatsky.utils.get_angle(nx,ny,D)
  """
  return flatsky.utils.get_angle(nx,ny,D)

def cutmap(ox,oy,cx,cy,omap):
  """
  Usage:
    - e.g., cmap = flatsky.utils.cutmap(ox,oy,cx,cy,omap)
  """
  return flatsky.utils.cutmap(ox,oy,cx,cy,omap)

