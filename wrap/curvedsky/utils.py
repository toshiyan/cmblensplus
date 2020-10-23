import libcurvedsky

def gauss1alm(lmax,Cl):
  """
  Generating alm as a random Gaussian field whose power spectrum is cl. The output alm is given by a 2D array.

  Args:
    :lmax (*int*): Maximum multipole of the output alm
    :Cl [*l*] (*double*): Angular power spectrum, with bounds (0:lmax)

  Returns:
    :alm [*l,m*] (*dcmplx*): Random Gaussian alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss1alm(lmax,Cl):
  """
  return libcurvedsky.utils.gauss1alm(lmax,Cl)

def gauss2alm(lmax,cl1,cl2,xl,flm=None):
  """
  Generating two alms as random Gaussian fields whose power spectra are cl1, cl2 and the cross spectrum is xl.

  Args:
    :lmax (*int*): Maximum multipole of the output alm
    :cl1 [*l*] (*double*): Angular power spectrum of the 1st alm, with bounds (0:lmax)
    :cl2 [*l*] (*double*): Angular power spectrum of the 2nd alm, with bounds (0:lmax)
    :xl [*l*] (*double*): Cross-power spectrum between alm1 and alm2, with bounds (0:lmax)

  Args(optional):
    :flm [*l,m*] (*dcmplx*): pre-computed Gaussian fields whose spectrum is cl1, with bounds (0:lmax,0:lmax), default to None

  Returns:
    :alm [*2,l,m*] (*dcmplx*): Random Gaussian alms, with bounds (2,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss2alm(lmax,cl1,cl2,xl,flm):
  """
  if flm is None:  flm = 0
  return libcurvedsky.utils.gauss2alm(lmax,cl1,cl2,xl,flm)

def gaussTEB(lmax,TT,EE,BB,TE):
  """
  Generating T/E/B alms as random Gaussian fields whose power spectra are TT, EE, BB and the cross spectrum is TE.

  Args:
    :lmax (*int*): Maximum multipole of the output alms
    :TT [*l*] (*double*): Angular power spectrum of temperature, with bounds (0:lmax)
    :EE [*l*] (*double*): Angular power spectrum of E mode, with bounds (0:lmax)
    :BB [*l*] (*double*): Angular power spectrum of B mode, with bounds (0:lmax)
    :TE [*l*] (*double*): TE cross-power spectrum, with bounds (0:lmax)

  Returns:
    :alm [*3,l,m*] (*dcmplx*): Random Gaussian T/E/B alms, with bounds (3,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gaussTEB(lmax,TT,EE,BB,TE):
  """
  return libcurvedsky.utils.gaussTEB(lmax,TT,EE,BB,TE)

def gauss3alm(lmax,cl):
  """
  Generating three alms as random Gaussian fields whose covariance is given by cl[i,j].

  Args:
    :lmax (*int*): Maximum multipole of the output alm
    :cl [*i,j,l*] (*double*): Covariance between the gaussian fields, with bounds (3,3,0:lmax)

  Returns:
    :alm [*3,l,m*] (*dcmplx*): Random Gaussian alms, with bounds (3,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss3alm(lmax,cl):
  """
  return libcurvedsky.utils.gauss3alm(lmax,cl)

def gauss4alm(lmax,cl):
  """
  Generating four alms as random Gaussian fields whose covariance is given by cl[i,j].

  Args:
    :lmax (*int*): Maximum multipole of the output alm
    :cl [*i,j,l*] (*double*): Covariance between the gaussian fields, with bounds (4,4,0:lmax)

  Returns:
    :alm [*4,l,m*] (*dcmplx*): Random Gaussian alms, with bounds (4,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss4alm(lmax,cl):
  """
  return libcurvedsky.utils.gauss4alm(lmax,cl)

def get_baseline(npix,nside_subpatch,QU):
  """
  Calculate baseline of each subpatch. The subpatches have the same size.
  Written by Ryo Nagata.

  Args:
    :npix (*int*): pixel number of the full map
    :nside_subpatch (*int*): Nside of sub patch
    :QU [*pix,2*] (*double*): Q/U maps, with bounds (0:npix-1,2)

  Returns:
    :blmap [*pix,2*] (*double*): baseline maps, with bounds (0:npix-1,2)

  Usage:
    :blmap = curvedsky.utils.get_baseline(npix,nside_subpatch,QU):
  """
  return libcurvedsky.utils.get_baseline(npix,nside_subpatch,QU)

def get_winmap(nside_large,nside_small,ipix_pix,apod):
  """
  Return apodization window for subpatch.
  Written by Ryo Nagata.

  Args:
    :nside_large (*int*): Nside of sub patch
    :nside_small (*int*): full Nside
    :ipix_pix (*int*): pixel index of full map
    :apod (*double*): apodization length

  Returns:
    :wind_out (*double*): aporization window at ipix_pix

  Usage:
    :win_out = curvedsky.utils.get_winmap(nside_large,nside_small,ipix_pix,apod):
  """
  return libcurvedsky.utils.get_winmap(nside_large,nside_small,ipix_pix,apod)

def get_apod_window(s,a):
  """
  A sine apodization window

  Args:
    :s (*double*): Distance from the center of the window
    :a (*double*): Apodization length, nothing (a=1) to all (a=0)

  Returns:
    :w (*double*): Aporization window

  Usage:
    :w = curvedsky.utils.get_apod_window(s,a):
  """
  return libcurvedsky.utils.get_apod_window(s,a)

def eb_separate(nside,lmax,W,Q,U):
  """
  E/B mode seperation based on the chi-field estimator. See e.g. Sec.III.2 of arXiv:1305.7441 for numerical implimentation.

  Args:
    :npix (*int*): Pixel number of the desired map
    :lmax (*int*): Maximum multipole used for the harmonic transform internally
    :W[*pix*] (*double*): Window function satisfying zero 1st and 2nd derivatives at the boundary, with bounds (0:npix-1)
    :Q/U[*pix*] (*double*): Input Q/U map already multiplied by W, with bounds (0:npix-1)

  Returns:
    :Elm/Blm[*l,m*] (*dcmplx*): Seperated E/B modes, with bounds (0:lmax,0:lmax)

  Usage:
    :Elm,Blm = curvedsky.utils.eb_separate(nside,lmax,W,Q,U):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.eb_separate(npix,lmax,W,Q,U)

def alm2cl(lmax,alm1,alm2=None):
  """
  From alm to angular power spectrum

  Args:
    :lmax (*int*): Maximum multipole of the input alm
    :alm1 [*l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)

  Args(optional):
    :alm2 [*l,m*] (*dcmplx*): 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1

  Returns:
    :cl [*l*] (*double*): Auto or cross angular power spectrum, with bounds (0:lmax)

  Usage:
    :cl = curvedsky.utils.alm2cl(lmax,alm1,alm2):
  """
  if alm2 is None:  alm2 = alm1
  return libcurvedsky.utils.alm2cl(lmax,alm1,alm2)

def alm2bcl(bn,lmax,alm1,spc='',alm2=None):
  """
  From alm to angular power spectrum with multipole binning

  Args:
    :bn (*int*): Number of multipole bins
    :lmax (*int*): maximum multipole of the input alm
    :alm1 [*l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)

  Args(optional):
    :alm2 [*l,m*] (*dcmplx*): 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1
    :spc (*str*): Specify bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing

  Returns:
    :cb [*bin*] (*double*): Auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)

  Usage:
    :cb = curvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,spc):
  """
  if alm2 is None:  alm2 = alm1
  return libcurvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,spc)

def alm2rho(lmax,alm1,alm2):
  """
  Compute correlation coefficients between two alms

  Args:
    :lmax (*int*): Maximum multipole of the input alm
    :alm1 [*l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
    :alm2 [*l,m*] (*dcmplx*): 2nd harmonic coefficient, with bounds (0:lmax,0:lmax)

  Returns:
    :rho [*l*] (*double*): Auto or cross angular power spectrum, with bounds (0:lmax)

  Usage:
    :rho = curvedsky.utils.alm2rho(lmax,alm1,alm2):
  """
  return libcurvedsky.utils.alm2rho(lmax,alm1,alm2)

def alm2cov(alm,n=0,lmax=0):
  """
  Compute correlation coefficients between two alms

  Args:
    :lmax (*int*): Maximum multipole of the input alms
    :alm [*n,l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)

  Returns:
    :cov [*n,n,l*] (*double*): Auto and cross angular power spectra between alm[*i*] and alm[*j*]

  Usage:
    :cov = curvedsky.utils.alm2cov(n,lmax,alm):
  """
  n = len(alm[:,0,0])
  lmax = len(alm[0,:,0]) - 1
  return libcurvedsky.utils.alm2cov(n,lmax,alm)

def apodize(nside,rmask,ascale,order=1,holeminsize=0):
  """
  Compute apodized window function. Partially use Healpix's process_mask code.

  Args:
    :nside (*int*): Nside of the input map
    :rmask[*pix*] (*double*): Input window function, with bounds (0:pix-1). Pixels at rmask=0 is considered as masked pixels.
    :ascale (*double*): Apodization length [*deg*] from the closest masked pixel

  Args(optional):
    :order (*int*): Pixel order, 1 for RING (default), otherwize NESTED
    :holeminsize (*double*): Minimum hole size [*arcmin*] (i.e., holes within this size is filled), default to 0

  Returns:
    :amask[*pix*] (*double*): Apodization window, with bounds (0:npix-1), using the same ordering as input

  Usage:
    :amask = curvedsky.utils.apodize(nside,rmask,ascale,order,holeminsize):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.apodize(npix,rmask,ascale,order,holeminsize)

def almxfl(lmax,mmax,alm,cl):
  """
  Calculate xlm[l,m] = alm[l,m] x cl[l]

  Args:
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:lmax)
    :cl [*l*] (*double*): 1D spectrum to be multiplied to alm, with bounds (0:lmax)

  Returns:
    :xlm [*l,m*] (*dcmplx*): Modified alm, with bounds (0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.utils.almxfl(lmax,mmax,alm,cl):
  """
  return libcurvedsky.utils.almxfl(lmax,mmax,alm,cl)

def hp_alm2map(nside,lmax,mmax,alm):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :nside (*int*): Nside of the input map
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)

  Returns:
    :map [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Usage:
    :map = curvedsky.utils.hp_alm2map(nside,lmax,mmax,alm):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.hp_alm2map(npix,lmax,mmax,alm)

def hp_alm2map_spin(nside,lmax,mmax,spin,elm,blm):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :nside (*int*): Nside of the input map
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :spin (*int*): Spin of the transform
    :elm [*l,m*] (*dcmplx*): Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
    :blm [*l,m*] (*dcmplx*): Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)

  Returns:
    :map0 [*pix*] (*double*): Real part of the transformed map (Q-like map), with bounds (0:npix-1)
    :map1 [*pix*] (*double*): Imaginary part of the transformed map (U-like map), with bounds (0:npix-1)

  Usage:
    :map0,map1 = curvedsky.utils.hp_alm2map_spin(nside,lmax,mmax,spin,elm,blm):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.hp_alm2map_spin(npix,lmax,mmax,spin,elm,blm)

def hp_map2alm(nside,lmax,mmax,map):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :nside (*int*): Nside of the input map
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :map [*pix*] (*double*): Input map, with bounds (0:npix-1)

  Returns:
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient obtained from the input map, with bounds (0:lmax,0:mmax)

  Usage:
    :alm = curvedsky.utils.hp_map2alm(nside,lmax,mmax,map):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.hp_map2alm(npix,lmax,mmax,map)

def hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1):
  """
  Spin Ylm transform of the map ( = map0 + i map1 ) to alm with the healpix (l,m) order. For example, if map0=Q, map1=U and spin=2, 
  the alm contains E-mode and B-mode. 

  Args:
    :nside (*int*): Nside of the input map
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :spin (*int*): Spin of the transform
    :map0 [*pix*] (*double*): Real part of the input map, with bounds (0:npix-1)
    :map1 [*pix*] (*double*): Imaginary part of the input map, with bounds (0:npix-1)

  Returns:
    :alm [*2,l,m*] (*dcmplx*): Parity-eve/odd harmonic coefficients obtained from the input map, with bounds (2,0:lmax,0:mmax)

  Usage:
    :alm = curvedsky.utils.hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.hp_map2alm_spin(npix,lmax,mmax,spin,map0,map1)

def map_mul_lfunc(nside,imap,lmax,lfunc):
  """
  Convert map to alm, multiply a function to alm and convert back again to map

  Args:
    :nside (*int*): Nside of input map
    :imap [*pix*] (*double*): Input map, with bounds (0:12*nside**2-1)
    :lmax (*int*): Maximum multipole of the input alm
    :lfunc [*l*] (*double*): 1D spectrum to be multiplied to alm, with bounds (0:lmax)

  Returns:
    :omap [*pix*] (*double*): Output map, with bounds (0:12*nside**2-1)

  Usage:
    :omap = curvedsky.utils.map_mul_lfunc(nside,imap,lmax,lfunc):
  """
  return libcurvedsky.utils.map_mul_lfunc(nside,imap,lmax,lfunc)

def mulwin(alm,win,nside=0,lmax=0,mmax=0):
  """
  Multiply window to a map obtained from alm

  Args:
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient to be multiplied at window, with bounds (0:lmax,0:mmax)
    :win [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Returns:
    :wlm [*l,m*] (*dcmplx*): Harmonic coefficient of the window-multiplied map, with bounds (0:lmax,0:mmax)

  Usage:
    :wlm = curvedsky.utils.mulwin(nside,lmax,mmax,alm,win):
  """
  npix = len(win)
  lmax = len(alm[:,0]) - 1
  mmax = len(alm[0,:]) - 1
  return libcurvedsky.utils.mulwin(npix,lmax,mmax,alm,win)

def mulwin_spin(elm,blm,win,nside=0,lmax=0,mmax=0,spin=2):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :elm [*l,m*] (*dcmplx*): Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
    :blm [*l,m*] (*dcmplx*): Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
    :win [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Args (Optional):
    :spin (*int*): Spin of the transform, default to 2

  Returns:
    :wlm [*2,l,m*] (*dcmplx*): Parity-eve/odd harmonic coefficients obtained from the window-multiplied map, with bounds (2,0:lmax,0:mmax)

  Usage:
    :wlm = curvedsky.utils.mulwin_spin(nside,lmax,mmax,spin,elm,blm,win):
  """
  npix = len(win)
  lmax = len(elm[:,0]) - 1
  mmax = len(elm[0,:]) - 1
  return libcurvedsky.utils.mulwin_spin(npix,lmax,mmax,spin,elm,blm,win)

def lm_healpy2healpix(almpy,lmax,lmpy=0):
  """
  Transform healpy alm to healpix alm

  Args:
    :lmax (*int*): Maximum multipole of the input/output alm satisfying 2 x lmpy = (lmax+1) x (lmax+2)
    :almpy[*index*] (*dcmplx*): Healpy alm, with bounds (0:lmpy-1)

  Returns:
    :almpix [*l,m*] (*dcmplx*): Healpix alm, with bounds (0:lmax,0:lmax)

  Usage:
    :almpix = curvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax):
  """
  lmpy = len(almpy)
  return libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax)

def cosin_healpix(nside):
  """
  Return cos(theta) as a function of the Healpix pixel index

  Args:
    :nside (*int*): Nside of the desired map

  Returns:
    :cosin [*pix*] (*double*): cosin(theta), with bounds (0:npix-1)

  Usage:
    :cosin = curvedsky.utils.cosin_healpix(nside):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.cosin_healpix(npix)

def load_defangle_takahashi(fname,npix,verbose=False):
  """
  Usage:
    :theta,phi = curvedsky.utils.load_defangle_takahashi(fname,npix,verbose):
  """
  return libcurvedsky.utils.load_defangle_takahashi(fname,npix,verbose)

def polcoord2angle(npix,theta,phi,verbose=False):
  """
  Usage:
    :angle = curvedsky.utils.polcoord2angle(npix,theta,phi,verbose):
  """
  return libcurvedsky.utils.polcoord2angle(npix,theta,phi,verbose)

def polcoord2angle_alm(nside,lmax,theta,phi,verbose=False):
  """
  Usage:
    :glm,clm = curvedsky.utils.polcoord2angle_alm(nside,lmax,theta,phi,verbose):
  """
  return libcurvedsky.utils.polcoord2angle_alm(nside,lmax,theta,phi,verbose)

def calc_mfs(bn,nu,lmax,walm,nside=0):
  """
  Compute 2D Minkowski functionals

  Args:
    :bn (*int*): Number of nu bins
    :nu [*bin*] (*double*): Nu bins, with bounds (bn)
    :lmax (*int*): Maximum multipole of the input walm
    :walm [*l,m*] (*dcmplx*): Alm with filtering, possibly divided by the map variance, with bounds (0:lmax,0:lmax)

  Args(optional): 
    :nside (*int*): Nside of the intermediate map, default to lmax

  Returns:
    :V [*bin,type*] (*double*): The three Minkowski functionals, V0, V1 and V2, at each nu bin, with bounds (bn,0:2)

  Usage:
    :V = curvedsky.utils.calc_mfs(bn,nu,lmax,walm,nside):
  """
  if nside == 0:  nside = lmax
  return libcurvedsky.utils.calc_mfs(bn,nu,lmax,walm,nside)

def mock_galaxy_takahashi(fname,zn,ngz,zi,a=2.0,b=1.0,zm=1.0,sz=0.0,zbias=0.0,b0=1.0,btype='sqrtz'):
  """
  Compute galaxy overdensity map from dark matter density map

  Args:
    :fname (*str*): Filename of density map
    :zn (*int*): Number of tomographic bins
    :ngz[*zn*] (*double*): Total number of galaxies at each z-bin
    :zi[*zn+1*] (*double*): Redshift intervals of the tomographic bins

  Args(optional): 
    :a (*double*): galaxy distribution shape parameter
    :b (*double*): galaxy distribution shape parameter
    :zm (*double*): mean redshift of galaxy distribution
    :b0 (*double*): constant galaxy bias at z=0

  Returns:
    :V [*bin,type*] (*double*): The three Minkowski functionals, V0, V1 and V2, at each nu bin, with bounds (bn,0:2)

  Usage:
    :gmap = curvedsky.utils.mock_galaxy_takahashi(fname,zn,ngz,zi,b0,btype,a,b,zm,sz,zbias):
  """
  return libcurvedsky.utils.mock_galaxy_takahashi(fname,zn,ngz,zi,b0,btype,a,b,zm,sz,zbias)

