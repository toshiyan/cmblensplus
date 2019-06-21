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

def gauss2alm(lmax,cl1,cl2,xl):
  """
  Generating two alms as random Gaussian fields whose power spectra are cl1, cl2 and the cross spectrum is xl.

  Args:
    :lmax (*int*): Maximum multipole of the output alm
    :cl1 [*l*] (*double*): Angular power spectrum of the 1st alm, with bounds (0:lmax)
    :cl2 [*l*] (*double*): Angular power spectrum of the 2nd alm, with bounds (0:lmax)
    :xl [*l*] (*double*): Cross-power spectrum between alm1 and alm2, with bounds (0:lmax)

  Returns:
    :alm [*2,l,m*] (*dcmplx*): Random Gaussian alms, with bounds (2,0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.utils.gauss2alm(lmax,cl1,cl2,xl):
  """
  return libcurvedsky.utils.gauss2alm(lmax,cl1,cl2,xl)

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
  Generating three alms as random Gaussian fields whose covariance is given by cl[*i,j*].

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
  Generating four alms as random Gaussian fields whose covariance is given by cl[*i,j*].

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

def eb_separate(npix,lmax,W,Q,U):
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
    :Elm,Blm = curvedsky.utils.eb_separate(npix,lmax,W,Q,U):
  """
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
  if alm2 is None: alm2= alm1
  return libcurvedsky.utils.alm2cl(lmax,alm1,alm2)

def alm2bcl(bn,lmax,alm1,alm2=None,spc= ''):
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
  if alm2 is None: alm2= alm1
  return libcurvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,spc)

def apodize(npix,rmask,ascale,order= 1,holeminsize= 0):
  """
  Compute apodized window function. Partially use Healpix's process_mask code.

  Args:
    :npix (*int*): Number of pixel
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
  return libcurvedsky.utils.apodize(npix,rmask,ascale,order,holeminsize)

def hp_alm2map(npix,lmax,mmax,alm):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :npix (*int*): Pixel number of the input map
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :alm [*l,m*] (*dcmplx*): Harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:lmax)

  Returns:
    :map [*pix*] (*double*): Transformed map, with bounds (0:npix-1)

  Usage:
    :map = curvedsky.utils.hp_alm2map(npix,lmax,mmax,alm):
  """
  return libcurvedsky.utils.hp_alm2map(npix,lmax,mmax,alm)

def hp_alm2map_spin(npix,lmax,mmax,spin,elm,blm):
  """
  Ylm transform of the map to alm with the healpix (l,m) order

  Args:
    :npix (*int*): Pixel number of the input map
    :lmax (*int*): Maximum multipole of the input alm
    :mmax (*int*): Maximum m of the input alm
    :spin (*int*): Spin of the transform
    :elm [*l,m*] (*dcmplx*): Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:lmax)
    :blm [*l,m*] (*dcmplx*): Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:lmax)

  Returns:
    :map0 [*pix*] (*double*): Real part of the transformed map (Q-like map), with bounds (0:npix-1)
    :map1 [*pix*] (*double*): Imaginary part of the transformed map (U-like map), with bounds (0:npix-1)

  Usage:
    :map0,map1 = curvedsky.utils.hp_alm2map_spin(npix,lmax,mmax,spin,elm,blm):
  """
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
    :alm1 [*l,m*] (*dcmplx*): Harmonic coefficient obtained from the input map, with bounds (0:lmax,0:mmax)

  Usage:
    :alm = curvedsky.utils.hp_map2alm(nside,lmax,mmax,map):
  """
  return libcurvedsky.utils.hp_map2alm(nside,lmax,mmax,map)

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
  return libcurvedsky.utils.hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1)

def lm_healpy2healpix(lmpy,almpy,lmax):
  """
  Transform healpy alm to healpix alm 

  Args:
    :lmpy (*int*): Length of healpy alm
    :lmax (*int*): Maximum multipole
    :almpy[*index*] (*dcmplx*): Healpy alm, with bounds (0:lmpy-1)

  Returns:
    :almpix [*l,m*] (*dcmplx*): Healpix alm, with bounds (0:lmax,0:lmax)

  Usage:
    :almpix = curvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax):
  """
  return libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax)

def cosin_healpix(npix,lmax):
  """
  Return cos(theta) as a function of the Healpix pixel index

  Args:
    :npix (*int*): Pixel number of the desired map
    :lmax (*int*): Maximum multipole

  Returns:
    :cosin [*pix*] (*double*): cosin(theta), with bounds (0:npix-1)

  Usage:
    :cosin = curvedsky.utils.cosin_healpix(npix,lmax):
  """
  return libcurvedsky.utils.cosin_healpix(npix,lmax)

def calc_mfs(bn,nu,lmax,walm,nside=None):
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
  if nside is None: nside= lmax
  return libcurvedsky.utils.calc_mfs(bn,nu,lmax,walm,nside)

def cinv(npix,lmax,cl,nij,alm,itern):
  """
 Computing inverse-variance filtered multipoles: C^-1d
 The power spectrum is assumed to be an isotropic Gaussian spectrum. 
 Inverse noise covariance is given in pixel space and uncorrelated (nij = sigma x delta_ij).

  Args:
    :npix (*int*): Number of pixel
    :lmax (*int*): Maximum multipole of alm
    :cl[*l*] (*double*): Angular power spectrum of alm, with bounds (0:lmax)
    :nij[*pix*] (*double*): Inverse of the noise variance at each pixel, with bounds (0:npix-1)
    :alm[*l,m*] (*dcmplx*): Input alm, with bouds (0:lmax,0:lmax)
    :itern (*int*): Number of interation

  Returns:
    :xlm[*l,m*] (*dcmplx*): C-inverse filtered multipoles, with bounds (0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.utils.cinv(npix,lmax,cl,nij,alm,itern):
  """
  return libcurvedsky.utils.cinv(npix,lmax,cl,nij,alm,itern)

def cg_algorithm(npix,lmax,clh,nij,b,itern):
  """
 Searching a solution of Ax=b with the Conjugate Gradient Algorithm iteratively
  Usage:
    :x = curvedsky.utils.cg_algorithm(npix,lmax,clh,nij,b,itern):
  """
  return libcurvedsky.utils.cg_algorithm(npix,lmax,clh,nij,b,itern)

def mat_multi(npix,lmax,clh,nij,x,mtype= ''):
  """
 multiplying matrix
  Usage:
    :v = curvedsky.utils.mat_multi(npix,lmax,clh,nij,x,mtype):
  """
  return libcurvedsky.utils.mat_multi(npix,lmax,clh,nij,x,mtype)

