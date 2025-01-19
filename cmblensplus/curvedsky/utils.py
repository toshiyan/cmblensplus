from cmblensplus import libcurvedsky
import numpy

def gauss1alm(cl,lmin=2):
  """
  Generating alm as a random Gaussian field whose power spectrum is cl. The output alm is given by a 2D array.

  Args:
    :cl [*l*] (*double*): Angular power spectrum, with bounds (0:lmax)

  Args(optional):
    :lmin (*int*): Minimum multipole of output alm (default: 2)

  Returns:
    :alm [*l,m*] (*dcmplx*): Random Gaussian alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss1alm(lmax,cl,lmin):
  """
  lmax = len(cl) - 1
  if cl[lmin] <= 0: raise SystemExit('ERROR in "gauss1alm": cl[lmin] should be positive')
  return libcurvedsky.utils.gauss1alm(lmax,cl,lmin)

def gauss2alm(cl1,cl2,xl,lmin=2,flm=None):
  """
  Generating two alms as random Gaussian fields whose power spectra are cl1, cl2 and the cross spectrum is xl.

  Args:
    :cl1 [*l*] (*double*): Angular power spectrum of the 1st alm, with bounds (0:lmax)
    :cl2 [*l*] (*double*): Angular power spectrum of the 2nd alm, with bounds (0:lmax)
    :xl [*l*] (*double*): Cross-power spectrum between alm1 and alm2, with bounds (0:lmax)

  Args(optional):
    :lmin (*int*): Minimum multipole of output alm (default: 2)
    :flm [*l,m*] (*dcmplx*): Constrained realiation of alms whose power spectrum is cl1

  Returns:
    :alm [*2,l,m*] (*dcmplx*): Random Gaussian alms, with bounds (2,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss2alm(lmax,cl1,cl2,xl,lmin,flm):
  """
  if flm is None: flm = numpy.zeros((lmax+1,lmax+1),dtype=numpy.complex128)
  lmax = len(cl1) - 1
  if numpy.any(xl[lmin:]**2/cl1[lmin:]/cl2[lmin:]>1): raise SystemExit('ERROR in "gauss2alm": correlation coefficient is larger than 1')
  return libcurvedsky.utils.gauss2alm(lmax,cl1,cl2,xl,lmin,flm)

def gaussTEB(TT,EE,BB,TE,lmin=2):
  """
  Generating T/E/B alms as random Gaussian fields whose power spectra are TT, EE, BB and the cross spectrum is TE.

  Args:
    :TT [*l*] (*double*): Angular power spectrum of temperature, with bounds (0:lmax)
    :EE [*l*] (*double*): Angular power spectrum of E mode, with bounds (0:lmax)
    :BB [*l*] (*double*): Angular power spectrum of B mode, with bounds (0:lmax)
    :TE [*l*] (*double*): TE cross-power spectrum, with bounds (0:lmax)

  Args(optional):
    :lmin (*int*): Minimum multipole of output alm (default: 2)

  Returns:
    :alm [*3,l,m*] (*dcmplx*): Random Gaussian T/E/B alms, with bounds (3,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gaussTEB(lmax,TT,EE,BB,TE,lmin):
  """
  lmax = len(TT) - 1
  if numpy.any(TT[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": TT[lmin:] should be positive')
  if numpy.any(EE[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": EE[lmin:] should be positive')
  if numpy.any(BB[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": BB[lmin:] should be positive')
  if numpy.any(TE[2:]**2/TT[2:]/EE[2:]>1): raise SystemExit('ERROR in "gaussTEB": correlation coefficient is larger than 1')
  return libcurvedsky.utils.gaussteb(lmax,TT,EE,BB,TE,lmin)

def gaussalm(cl,n=None,lmax=None,ilm=None):
  """
  Generating alms as random Gaussian fields whose covariance is given by cl[i,j].

  Args:
    :cl [*i,j,l*] (*double*): Covariance between the gaussian fields, with bounds (n,n,0:lmax)

  Args(optional):
    :ilm [*l,m*] (*dcmplx*): Input alm for the cl[*0,0*] element (default to None). The other alms are generated to be correlated with ilm.

  Returns:
    :alm [*i,l,m*] (*dcmplx*): Random Gaussian alms, with bounds (n,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gaussalm(n,lmax,cl,ilm):
  """
  if n is None:    n    = len(cl[:,0,0])
  if lmax is None: lmax = len(cl[0,0,:]) - 1
  if ilm is None:  ilm  = [[0 for x in range(lmax+1)] for y in range(lmax+1)] 
  return libcurvedsky.utils.gaussalm(n,lmax,cl,ilm)

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

def eb_separate(lmax,W,Q,U):
  """
  E/B mode seperation based on the chi-field estimator. See e.g. Sec.III.2 of arXiv:1305.7441 for numerical implimentation.

  Args:
    :lmax (*int*): Maximum multipole used for the harmonic transform internally
    :W[*pix*] (*double*): Window function satisfying zero 1st and 2nd derivatives at the boundary, with bounds (0:npix-1)
    :Q/U[*pix*] (*double*): Input Q/U map already multiplied by W, with bounds (0:npix-1)

  Returns:
    :Elm/Blm[*l,m*] (*dcmplx*): Seperated E/B modes, with bounds (0:lmax,0:lmax)

  Usage:
    :Elm,Blm = curvedsky.utils.eb_separate(npix,lmax,W,Q,U):
  """
  npix = len(W)
  return libcurvedsky.utils.eb_separate(npix,lmax,W,Q,U)

def alm2cl(alm1,alm2=None):
  """
  From alm to angular power spectrum

  Args:
    :alm1 [*l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)

  Args(optional):
    :alm2 [*l,m*] (*dcmplx*): 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1

  Returns:
    :cl [*l*] (*double*): Auto or cross angular power spectrum, with bounds (0:lmax)

  Usage:
    :cl = curvedsky.utils.alm2cl(lmax,alm1,alm2):
  """
  lmax = len(alm1[:,0]) - 1
  if alm2 is None:  alm2 = alm1
  return libcurvedsky.utils.alm2cl(lmax,alm1,alm2)

def alm2bcl(bn,alm1,alm2=None,spc=''):
  """
  From alm to angular power spectrum with multipole binning

  Args:
    :bn (*int*): Number of multipole bins
    :alm1 [*l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)

  Args(optional):
    :alm2 [*l,m*] (*dcmplx*): 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1
    :spc (*str*): Specify bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing

  Returns:
    :cb [*bin*] (*double*): Auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)

  Usage:
    :cb = curvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,spc):
  """
  lmax = len(alm1[:,0]) - 1
  if alm2 is None:  alm2 = alm1
  return libcurvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,spc)

def alm2rho(alm1,alm2):
  """
  Compute correlation coefficients between two alms

  Args:
    :alm1 [*l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
    :alm2 [*l,m*] (*dcmplx*): 2nd harmonic coefficient, with bounds (0:lmax,0:lmax)

  Returns:
    :rho [*l*] (*double*): Auto or cross angular power spectrum, with bounds (0:lmax)

  Usage:
    :rho = curvedsky.utils.alm2rho(lmax,alm1,alm2):
  """
  lmax = len(alm1[:,0]) - 1
  return libcurvedsky.utils.alm2rho(lmax,alm1,alm2)

def alm2cov(alm):
  """
  Compute correlation coefficients between two alms

  Args:
    :alm [*n,l,m*] (*dcmplx*): 1st harmonic coefficient, with bounds (0:lmax,0:lmax)

  Returns:
    :cov [*n,n,l*] (*double*): Auto and cross angular power spectra between alm[*i*] and alm[*j*]

  Usage:
    :cov = curvedsky.utils.alm2cov(n,lmax,alm):
  """
  n = len(alm[:,0,0])
  lmax = len(alm[0,:,0]) - 1
  return libcurvedsky.utils.alm2cov(n,lmax,alm)

def apodize(rmask,ascale,order=1,holeminsize=0):
  """
  Compute apodized window function. Partially use Healpix's process_mask code.

  Args:
    :rmask[*pix*] (*double*): Input window function, with bounds (0:pix-1). Pixels at rmask=0 is considered as masked pixels.
    :ascale (*double*): Apodization length [*deg*] from the closest masked pixel

  Args(optional):
    :order (*int*): Pixel order, 1 for RING (default), otherwize NESTED
    :holeminsize (*double*): Minimum hole size [*arcmin*] (i.e., holes within this size is filled), default to 0

  Returns:
    :amask[*pix*] (*double*): Apodization window, with bounds (0:npix-1), using the same ordering as input

  Usage:
    :amask = curvedsky.utils.apodize(npix,rmask,ascale,order,holeminsize):
  """
  npix = len(rmask)
  return libcurvedsky.utils.apodize(npix,rmask,ascale,order,holeminsize)

def hp_alm2map(nside,alm):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :nside (*int*): Nside of the input map
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)

  Returns:
    :map [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Usage:
    :map = curvedsky.utils.hp_alm2map(nside,lmax,mmax,alm):
  """
  lmax = len(alm[:,0]) - 1
  mmax = len(alm[0,:]) - 1
  npix = 12*nside**2
  return libcurvedsky.utils.hp_alm2map(npix,lmax,mmax,alm)

def hp_alm2map_spin(nside,spin,elm,blm):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :nside (*int*): Nside of the input map
    :spin (*int*): Spin of the transform
    :elm [*l,m*] (*dcmplx*): Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
    :blm [*l,m*] (*dcmplx*): Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)

  Returns:
    :map0 [*pix*] (*double*): Real part of the transformed map (Q-like map), with bounds (0:npix-1)
    :map1 [*pix*] (*double*): Imaginary part of the transformed map (U-like map), with bounds (0:npix-1)

  Usage:
    :map0,map1 = curvedsky.utils.hp_alm2map_spin(nside,lmax,mmax,spin,elm,blm):
  """
  lmax = len(elm[:,0]) - 1
  mmax = len(elm[0,:]) - 1
  npix = 12*nside**2
  return libcurvedsky.utils.hp_alm2map_spin(npix,lmax,mmax,spin,elm,blm)

def hp_map2alm(lmax,mmax,map):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :map [*pix*] (*double*): Input map, with bounds (0:npix-1)

  Returns:
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient obtained from the input map, with bounds (0:lmax,0:mmax)

  Usage:
    :alm = curvedsky.utils.hp_map2alm(npix,lmax,mmax,map):
  """
  npix = len(map)
  return libcurvedsky.utils.hp_map2alm(npix,lmax,mmax,map)

def hp_map2alm_spin(lmax,mmax,spin,map0,map1):
  """
  Spin Ylm transform of the map ( = map0 + i map1 ) to alm with the healpix (l,m) order. For example, if map0=Q, map1=U and spin=2, 
  the alm contains E-mode and B-mode. 

  Args:
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :spin (*int*): Spin of the transform
    :map0 [*pix*] (*double*): Real part of the input map, with bounds (0:npix-1)
    :map1 [*pix*] (*double*): Imaginary part of the input map, with bounds (0:npix-1)

  Returns:
    :alm [*2,l,m*] (*dcmplx*): Parity-eve/odd harmonic coefficients obtained from the input map, with bounds (2,0:lmax,0:mmax)

  Usage:
    :alm = curvedsky.utils.hp_map2alm_spin(npix,lmax,mmax,spin,map0,map1):
  """
  npix = len(map0)
  return libcurvedsky.utils.hp_map2alm_spin(npix,lmax,mmax,spin,map0,map1)

def map_mul_lfunc(imap,lfunc):
  """
  Convert map to alm, multiply a function to alm and convert back again to map

  Args:
    :imap [*pix*] (*double*): Input map, with bounds (0:12*nside**2-1)
    :lfunc [*l*] (*double*): 1D spectrum to be multiplied to alm, with bounds (0:lmax)

  Returns:
    :omap [*pix*] (*double*): Output map, with bounds (0:12*nside**2-1)

  Usage:
    :omap = curvedsky.utils.map_mul_lfunc(npix,imap,lmax,lfunc):
  """
  lmax = len(lfunc) - 1
  npix = len(imap)
  return libcurvedsky.utils.map_mul_lfunc(npix,imap,lmax,lfunc)

def mulwin(alm,win):
  """
  Multiply window to a map obtained from alm

  Args:
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient to be multiplied at window, with bounds (0:lmax,0:mmax)
    :win [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Returns:
    :wlm [*l,m*] (*dcmplx*): Harmonic coefficient of the window-multiplied map, with bounds (0:lmax,0:mmax)

  Usage:
    :wlm = curvedsky.utils.mulwin(npix,lmax,mmax,alm,win):
  """
  npix = len(win)
  lmax = len(alm[:,0]) - 1
  mmax = len(alm[0,:]) - 1
  return libcurvedsky.utils.mulwin(npix,lmax,mmax,alm,win)

def mulwin_spin(elm,blm,win,spin=2):
  """
  Multiply window to a map obtained from alm

  Args:
    :elm [*l,m*] (*dcmplx*): Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
    :blm [*l,m*] (*dcmplx*): Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
    :win [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Args (Optional):
    :spin (*int*): Spin of the transform, default to 2

  Returns:
    :wlm [*2,l,m*] (*dcmplx*): Parity-eve/odd harmonic coefficients obtained from the window-multiplied map, with bounds (2,0:lmax,0:mmax)

  Usage:
    :wlm = curvedsky.utils.mulwin_spin(npix,lmax,mmax,spin,elm,blm,win):
  """
  npix = len(win)
  lmax = len(elm[:,0]) - 1
  mmax = len(elm[0,:]) - 1
  return libcurvedsky.utils.mulwin_spin(npix,lmax,mmax,spin,elm,blm,win)

def lm_healpy2healpix(almpy):
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
  lmax = int((-3.+numpy.sqrt(1.+8.*lmpy))/2.)
  return libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax)

def cosin_healpix(nside):
  """
  Return cos(theta) as a function of the Healpix pixel index

  Args:
    :nside (*int*): Nside of the desired map

  Returns:
    :cosin[*pix*] (*double*): cosin(theta), with bounds (0:npix-1)

  Usage:
    :cosin = curvedsky.utils.cosin_healpix(nside):
  """
  npix = 12*nside**2
  return libcurvedsky.utils.cosin_healpix(npix)

def load_defangle_takahashi(fname,npix,verbose=False):
  """
  Read theta and phi coordinates at source plane obtained by Takahashi et al. (2017)

  Args:
    :fname (*str*): file name
    :npix (*int*): Number of pixels of theta and phi data

  Args (optional):
    :verbose (*bool*): output messages, default to False

  Returns:
    :theta[*pix*] (*double*): theta, with bounds (0:npix-1)
    :phi[*pix*] (*double*): phi, with bounds (0:npix-1)

  Usage:
    :theta,phi = curvedsky.utils.load_defangle_takahashi(fname,npix,verbose):
  """
  return libcurvedsky.utils.load_defangle_takahashi(fname,npix,verbose)

def polcoord2angle(npix,theta,phi,verbose=False):
  """
  Converting theta and phi coordinates at source plane to deflection angle.
  The algorithm is provided by Takashi Hamana and Ryuichi Takahashi.

  Args:
    :npix (*int*): Number of pixels of theta and phi data
    :theta[*pix*] (*double*): theta, with bounds (0:npix-1)
    :phi[*pix*] (*double*): phi, with bounds (0:npix-1)

  Args (optional):
    :verbose (*bool*): output messages, default to False

  Returns:
    :angle[*pix,2*] (*double*): deflection angle vector containing two components, with bounds (0:npix-1,1:2)

  Usage:
    :angle = curvedsky.utils.polcoord2angle(npix,theta,phi,verbose):
  """
  return libcurvedsky.utils.polcoord2angle(npix,theta,phi,verbose)

def polcoord2angle_alm(nside,lmax,theta,phi,verbose=False):
  """
  Converting theta and phi coordinates at source plane to deflection angle.

  Args:
    :npix (*int*): Number of pixels of theta and phi data
    :lmax (*int*): Maximum multipole of alms for gradient and curl modes
    :theta[*pix*] (*double*): theta, with bounds (0:npix-1)
    :phi[*pix*] (*double*): phi, with bounds (0:npix-1)

  Args (optional):
    :verbose (*bool*): output messages, default to False

  Returns:
    :glm[*l,m*] (*dcmplx*): gradient mode, with bounds (0:lmax,0:lmax)
    :clm[*l,m*] (*dcmplx*): curl mode, with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.utils.polcoord2angle_alm(nside,lmax,theta,phi,verbose):
  """
  return libcurvedsky.utils.polcoord2angle_alm(nside,lmax,theta,phi,verbose)

def calc_mfs(bn,nu,lmax,walm,nside=0):
  """
  Compute 2D Minkowski functionals from a given alm

  Args:
    :bn (*int*): Number of nu bins
    :nu [*bin*] (*double*): Nu bins, with bounds (bn)
    :lmax (*int*): Maximum multipole of the input walm
    :walm [*l,m*] (*dcmplx*): Alm with filtering, possibly divided by the map variance, with bounds (0:lmax,0:lmax)

  Args(optional): 
    :nside (*int*): Nside of the intermediate map, default to the closest power of 2 to 3xlmax

  Returns:
    :V [*bin,type*] (*double*): The three Minkowski functionals, V0, V1 and V2, at each nu bin, with bounds (bn,0:2)

  Usage:
    :V = curvedsky.utils.calc_mfs(bn,nu,lmax,walm,nside):
  """
  return libcurvedsky.utils.calc_mfs(bn,nu,lmax,walm,nside)

def mock_galaxy_takahashi(fname,zn,ngz,zi,b0=1.0,btype='sqrtz',a=2.0,b=1.0,zm=1.0,sz=0.0,zbias=0.0):
  """
  Compute galaxy overdensity map from dark matter density map of Takahashi et al. (2017). 
  The galaxy z distribution is assumed to have the functional form given by Eq.(7) of https://arxiv.org/abs/1810.03346

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
    :gmap [*pix,zbin*] (*double*): The galaxy number density map at each zbin

  Usage:
    :gmap = curvedsky.utils.mock_galaxy_takahashi(fname,zn,ngz,zi,b0,btype,a,b,zm,sz,zbias):
  """
  return libcurvedsky.utils.mock_galaxy_takahashi(fname,zn,ngz,zi,b0,btype,a,b,zm,sz,zbias)

