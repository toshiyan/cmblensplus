import curvedsky

def gauss1alm(lmax,cl):
  """
  Generating alm as a random Gaussian field whose power spectrum is cl.
  The output alm is given by a 2D array.

  Args:
    - lmax (int)        : maximum multipole of the output alm
    - Cl[l] (double)    : angular power spectrum, with bounds (1:lmax)

  Returns:
    - alm[l,m] (dcmplx) : random Gaussian alm, with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., alm = curvedsky.utils.gauss1alm(lmax,cl)
  """
  return curvedsky.utils.gauss1alm(lmax,cl)

def gauss2alm(lmax,cl1,cl2,xl):
  """
  Generating two alms as random Gaussian fields whose power spectra are cl1, cl2 and the cross spectrum is xl.

  Args:
    - lmax (int)      : maximum multipole of the output alm
    - cl1[l] (double) : angular power spectrum of the 1st alm, with bounds (1:lmax)
    - cl2[l] (double) : angular power spectrum of the 2nd alm, with bounds (1:lmax)
    - xl[l] (double)  : cross-power spectrum between alm1 and alm2, with bounds (1:lmax)

  Returns:
    - alm[2,l,m] (dcmplx): random Gaussian alms, with bounds (2,0:lmax,0:lmax)

  Usage:
    - e.g., alm = curvedsky.utils.gauss2alm(lmax,cl1,cl2,xl)
  """
  return curvedsky.utils.gauss2alm(lmax,cl1,cl2,xl)

def gaussTEB(lmax,TT,EE,BB,TE):
  """
  Generating T/E/B alms as random Gaussian fields whose power spectra are TT, EE, BB and the cross spectrum is TE.

  Args:
    - lmax (int)     : maximum multipole of the output alms
    - TT[l] (double) : angular power spectrum of temperature, with bounds (1:lmax)
    - EE[l] (double) : angular power spectrum of E mode, with bounds (1:lmax)
    - BB[l] (double) : angular power spectrum of B mode, with bounds (1:lmax)
    - TE[l] (double) : TE cross-power spectrum, with bounds (1:lmax)

  Returns:
    - alm[3,l,m] (dcmplx): random Gaussian T/E/B alms, with bounds (3,0:lmax,0:lmax)

  Usage:
    - e.g., alm = curvedsky.utils.gaussTEB(lmax,TT,EE,BB,TE)
  """
  return curvedsky.utils.gaussTEB(lmax,TT,EE,BB,TE)

def gauss3alm(lmax,cl):
  """
  Generating three alms as random Gaussian fields whose covariance is given by cl[i,j].

  Args:
    - lmax (int)         : maximum multipole of the output alm
    - cl[i,j,l] (double) : covariance between the gaussian fields, with bounds (3,3,1:lmax)

  Returns:
    - alm[3,l,m] (dcmplx): random Gaussian alms, with bounds (3,0:lmax,0:lmax)

  Usage:
    - e.g., alm = curvedsky.utils.gauss3alm(lmax,cl)
  """
  return curvedsky.utils.gauss3alm(lmax,cl)

def gauss4alm(lmax,cl):
  """
  Generating four alms as random Gaussian fields whose covariance is given by cl[i,j].

  Args:
    - lmax (int)         : maximum multipole of the output alm
    - cl[i,j,l] (double) : covariance between the gaussian fields, with bounds (4,4,1:lmax)

  Returns:
    - alm[4,l,m] (dcmplx): random Gaussian alms, with bounds (4,0:lmax,0:lmax)

  Usage:
    - e.g., alm = curvedsky.utils.gauss4alm(lmax,cl)
  """
  return curvedsky.utils.gauss4alm(lmax,cl)

def get_baseline(npix,nside_subpatch,QU):
  """
  Calculate baseline of each subpatch. The subpatches have the same size.
  Written by Ryo Nagata.

  Args:
    - npix (int)           : pixel number of the full map
    - nside_subpatch (int) : Nside of sub patch
    - QU[pix,2] (double)   : Q/U maps, with bounds (0:npix-1,2)

  Returns:
    - blmap[pix,2] (double): baseline maps, with bounds (0:npix-1,2)

  Usage:
    - e.g., blmap = curvedsky.utils.get_baseline(npix,nside_subpatch,QU)
  """
  return curvedsky.utils.get_baseline(npix,nside_subpatch,QU)

def get_winmap(nside_large,nside_small,ipix_pix,apod):
  """
  Return apodization window for subpatch.
  Written by Ryo Nagata.

  Args:
    - nside_large (int) : Nside of sub patch
    - nside_small (int) : full Nside
    - ipix_pix (int)    : pixel index of full map
    - apod (double)     : apodization length

  Returns:
    - wind_out (double) : aporization window at ipix_pix 

  Usage:
    - e.g., win_out = curvedsky.utils.get_winmap(nside_large,nside_small,ipix_pix,apod)
  """
  return curvedsky.utils.get_winmap(nside_large,nside_small,ipix_pix,apod)

def get_apod_window(s,a):
  """
  A sine apodization window

  Args:
    - s (double) : distance from the center of the window
    - a (double) : apodization length, nothing (a=1) to all (a=0)

  Returns:
    - w (double) : aporization window

  Usage:
    - e.g., w = curvedsky.utils.get_apod_window(s,a)
  """
  return curvedsky.utils.get_apod_window(s,a)

def cosin_healpix(npix,lmax) :
  """
  cos(theta) as a function of the Healpix pixel index

  Args:
    - npix (int) : pixel number of the desired map
    - lmax (int) : maximum multipole

  Returns:
    - cosin[pix] (double) : cos(theta), with bounds (0:npix-1)

  Usage:
    - e.g., cosin = curvedsky.utils.cosin_healpix(npix,lmax) 
  """
  return curvedsky.utils.cosin_healpix(npix,lmax) 

def eb_separate(npix,lmax,W,QUin):
  """
  E/B mode seperation based on the chi-field estimator.

  Args:
    - npix (int)           : pixel number of the desired map
    - lmax (int)           : maximum multipole used for the harmonic transform internally
    - W[pix,2] (double)    : window function satisfying zero 1st and 2nd derivatives at the boundary, with bounds (0:npix-1,2)
    - QUin[pix,2] (double) : input QU map, with bounds (0:npix-1,2)

  Returns:
    - QUou[pix,2] (double) : E/B separated QU map, with bounds (0:npix-1,2)

  Usage:
    - e.g., QUou = curvedsky.utils.eb_separate(npix,lmax,W,QUin)
  """
  return curvedsky.utils.eb_separate(npix,lmax,W,QUin)

def alm2cl(lmax,alm1,alm2=0,norm=0):
  """
  Usage:
    - e.g., cl = curvedsky.utils.alm2cl(lmax,alm1,alm2,norm)
  """
  if alm2==0: alm2= alm1
  if norm==0: norm= 1
  return curvedsky.utils.alm2cl(lmax,alm1,alm2,norm)

def alm2bcl(bn,lmax,alm1,alm2=0,oL=0,norm=0):
  """
  Usage:
    - e.g., cb = curvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,oL,norm)
  """
  if alm2==0: alm2= alm1
  if oL==0: oL= 0
  if norm==0: norm= 1
  return curvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,oL,norm)

def apodize(npix,rmask,ascale,order=0,holeminsize=0):
  """
  Usage:
    - e.g., amask = curvedsky.utils.apodize(npix,rmask,ascale,order,holeminsize)
  """
  if order==0: order= 1
  if holeminsize==0: holeminsize= 0
  return curvedsky.utils.apodize(npix,rmask,ascale,order,holeminsize)

def hp_map2alm(nside,lmax,mmax,map):
  """
  Usage:
    - e.g., alm = curvedsky.utils.hp_map2alm(nside,lmax,mmax,map)
  """
  return curvedsky.utils.hp_map2alm(nside,lmax,mmax,map)

def hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1):
  """
  Usage:
    - e.g., alm = curvedsky.utils.hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1)
  """
  return curvedsky.utils.hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1)

