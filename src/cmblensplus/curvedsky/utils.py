from ._core import sht_ducc0 as sht
from . import libcurvedsky
import numpy as np


def gauss1alm(cl, lmin=2):
    r"""
    Generate :math:`a_{lm}` as a random Gaussian field whose power spectrum is ``cl``.

    The output alm is given by a 2D array.

    Parameters
    ----------
    cl : array_like of float, shape (lmax + 1,)
        Angular power spectrum, with bounds ``0:lmax``.
    lmin : int, optional
        Minimum multipole of output alm. Default is 2.

    Returns
    -------
    alm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Random Gaussian alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    
    lmax = len(cl) - 1
    if cl[lmin] <= 0:
        raise SystemExit('ERROR in "gauss1alm": cl[lmin] should be positive')
    return libcurvedsky.utils.gauss1alm(lmax, cl, lmin)
    

def gauss2alm(cl1,cl2,xl,lmin=2,flm=None):
    r"""
    Generate two correlated alms as random Gaussian fields.
  
    The auto spectra of the two fields are ``cl1`` and ``cl2``, and their cross spectrum is ``xl``.
  
    Parameters
    ----------
    cl1 : array_like of float, shape (lmax + 1,)
        Angular power spectrum of the first alm, with bounds ``0:lmax``.
    cl2 : array_like of float, shape (lmax + 1,)
        Angular power spectrum of the second alm, with bounds ``0:lmax``.
    xl : array_like of float, shape (lmax + 1,)
        Cross-power spectrum between the two alms, with bounds ``0:lmax``.
    lmin : int, optional
        Minimum multipole of the output alms. Default is 2.
    flm : array_like of complex, optional
        Constrained realization of an alm whose power spectrum is ``cl1``.
  
    Returns
    -------
    alm : ndarray of complex, shape (2, lmax + 1, lmax + 1)
        Random Gaussian alms, with bounds ``(2, 0:lmax, 0:lmax)``.
    """

    if flm is None: 
        flm = np.zeros((lmax+1,lmax+1),dtype=np.complex128)
    lmax = len(cl1) - 1
    
    if np.any(xl[lmin:]**2/cl1[lmin:]/cl2[lmin:]>1): 
        raise SystemExit('ERROR in "gauss2alm": correlation coefficient is larger than 1')
    
    return libcurvedsky.utils.gauss2alm(lmax,cl1,cl2,xl,lmin,flm)

    
def gaussTEB(TT,EE,BB,TE,lmin=2):
    r"""
    Generate T/E/B alms as correlated random Gaussian fields.
  
    The auto spectra are ``TT``, ``EE``, and ``BB``.  The T-E cross spectrum is ``TE``.
  
    Parameters
    ----------
    TT : array_like of float, shape (lmax + 1,)
        Temperature angular power spectrum, with bounds ``0:lmax``.
    EE : array_like of float, shape (lmax + 1,)
        E-mode angular power spectrum, with bounds ``0:lmax``.
    BB : array_like of float, shape (lmax + 1,)
        B-mode angular power spectrum, with bounds ``0:lmax``.
    TE : array_like of float, shape (lmax + 1,)
        T-E cross-power spectrum, with bounds ``0:lmax``.
    lmin : int, optional
        Minimum multipole of the output alms. Default is 2.
  
    Returns
    -------
    alm : ndarray of complex, shape (3, lmax + 1, lmax + 1)
        Random Gaussian T/E/B alms, with bounds ``(3, 0:lmax, 0:lmax)``.
    """

    lmax = len(TT) - 1
    if np.any(TT[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": TT[lmin:] should be positive')
    if np.any(EE[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": EE[lmin:] should be positive')
    if np.any(BB[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": BB[lmin:] should be positive')
    if np.any(TE[2:]**2/TT[2:]/EE[2:]>1): raise SystemExit('ERROR in "gaussTEB": correlation coefficient is larger than 1')
    return libcurvedsky.utils.gaussteb(lmax,TT,EE,BB,TE,lmin)


def gaussalm(cl,n=None,lmax=None,ilm=None):
    r"""
    Generate correlated alms as random Gaussian fields.
  
    The covariance between fields is specified by ``cl[i, j]``.
  
    Parameters
    ----------
    cl : array_like of float, shape (n, n, lmax + 1)
        Covariance spectra between the Gaussian fields, with bounds ``(n, n, 0:lmax)``.
    n : int, optional
        Number of fields. If omitted, it is inferred from ``cl``.
    lmax : int, optional
        Maximum multipole. If omitted, it is inferred from ``cl``.
    ilm : array_like of complex, optional
        Input fixed alm for the ``cl[0, 0]`` element. The other alms are generated to be correlated with ``ilm``.
  
    Returns
    -------
    alm : ndarray of complex, shape (n, lmax + 1, lmax + 1)
        Random Gaussian alms, with bounds ``(n, 0:lmax, 0:lmax)``.
    """
    
    if n is None:    n    = len(cl[:,0,0])
    if lmax is None: lmax = len(cl[0,0,:]) - 1
    if ilm is None:  ilm  = [[0 for x in range(lmax+1)] for y in range(lmax+1)] 
    return libcurvedsky.utils.gaussalm(n,lmax,cl,ilm)

    
def subpatch_mask(nside_full,nside_sub,ipix_sub,ascale):
    r"""
    Calculate a mask for a patched map.
    
    Parameters
    ----------
    nside_full : int
        Nside of the full map.
    nside_sub : int
        Nside of the patch map.
    ipix_sub : int
        Healpix pixel index of the patch map.
    ascale : float
        Apodization scale, from no apodization at 1 to full apodization at 0.
    
    Returns
    -------
    mask : ndarray of float, shape (npix, 2)
        Mask map, with bounds ``(0:npix-1, 2)``.
    """
    
    return lib_utils.subpatch_mask(nside_full,nside_sub,ipix_sub,ascale)


def get_baseline(npix,nside_subpatch,QU):
    r"""
    Calculate the baseline of each subpatch.
    
    The subpatches are assumed to have the same size.
    
    Parameters
    ----------
    npix : int
        Number of pixels in the full map.
    nside_subpatch : int
        Nside of each subpatch.
    QU : array_like of float, shape (npix, 2)
        Q/U maps, with bounds ``(0:npix-1, 2)``.
    
    Returns
    -------
    blmap : ndarray of float, shape (npix, 2)
        Baseline maps, with bounds ``(0:npix-1, 2)``.
    """
    
    return lib_utils.get_baseline(npix,nside_subpatch,QU)

    
def get_winmap(nside_large,nside_small,ipix_pix,apod):
    r"""
    Return the apodization window value for a subpatch.
    
    Parameters
    ----------
    nside_large : int
        Nside of the subpatch.
    nside_small : int
        Full-map Nside.
    ipix_pix : int
        Pixel index of the full map.
    apod : float
        Apodization length.
    
    Returns
    -------
    wind_out : float
        Apodization window value at ``ipix_pix``.
    """
    
    return lib_utils.get_winmap(nside_large,nside_small,ipix_pix,apod)

    
def get_apod_window(s,a):
    r"""
    Evaluate a sine apodization window.
    
    Parameters
    ----------
    s : float
        Distance from the center of the window.
    a : float
        Apodization length, from no apodization at 1 to full apodization at 0.
    
    Returns
    -------
    w : float
        Apodization window value.
    """
    
    return lib_utils.default_apod_window(s,a)

    
def eb_separate(lmax,W,Q,U):
    r"""
    Perform E/B-mode separation using the chi-field estimator.
    
    See Sec. III.2 of arXiv:1305.7441 for the numerical implementation.
    
    Parameters
    ----------
    lmax : int
        Maximum multipole used internally for the harmonic transform.
    W : array_like of float, shape (npix,)
        Window function satisfying zero first and second derivatives at the
        boundary, with bounds ``0:npix-1``.
    Q : array_like of float, shape (npix,)
        Input Q map already multiplied by ``W``, with bounds ``0:npix-1``.
    U : array_like of float, shape (npix,)
        Input U map already multiplied by ``W``, with bounds ``0:npix-1``.
    
    Returns
    -------
    Elm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Separated E modes, with bounds ``(0:lmax, 0:lmax)``.
    Blm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Separated B modes, with bounds ``(0:lmax, 0:lmax)``.
    """
    
    return sht.eb_separate(lmax,W,Q,U)

    
def alm2cl(alm1,alm2=None):
    r"""
    Compute an angular power spectrum from alm coefficients.
  
    Parameters
    ----------
    alm1 : array_like of complex, shape (lmax + 1, lmax + 1)
        First harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``.
    alm2 : array_like of complex, optional
        Second harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``. If omitted, ``alm1`` is used to compute the auto spectrum.
  
    Returns
    -------
    cl : ndarray of float, shape (lmax + 1,)
        Auto or cross angular power spectrum, with bounds ``0:lmax``.
    """

    lmax = len(alm1[:,0]) - 1
    if alm2 is None:  alm2 = alm1
    
    return libcurvedsky.utils.alm2cl(lmax,alm1,alm2)

    
def alm2bcl(bn,alm1,alm2=None,spc=''):
    r"""
    Compute a binned angular power spectrum from alm coefficients.
  
    Parameters
    ----------
    bn : int
        Number of multipole bins.
    alm1 : array_like of complex, shape (lmax + 1, lmax + 1)
        First harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``.
    alm2 : array_like of complex, optional
        Second harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``. If omitted, ``alm1`` is used.
    spc : str, optional
        Bin spacing. Use ``''`` for linear spacing, ``'log'`` for logarithmic spacing, 
        ``'log10'`` for base-10 logarithmic spacing, ``'p2'`` for powers of 2, 
        or ``'p3'`` for powers of 3.
  
    Returns
    -------
    cb : ndarray of float, shape (bn,)
        Auto or cross angular power spectrum with multipole binning, with bounds ``0:bn-1``.
    """
    
    lmax = len(alm1[:,0]) - 1
    if alm2 is None:  alm2 = alm1
    
    return libcurvedsky.utils.alm2bcl(bn,lmax,alm1,alm2,spc)


def alm2rho(alm1,alm2):
    r"""
    Compute correlation coefficients between two alms.
  
    Parameters
    ----------
    alm1 : array_like of complex, shape (lmax + 1, lmax + 1)
        First harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``.
    alm2 : array_like of complex, shape (lmax + 1, lmax + 1)
        Second harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``.
  
    Returns
    -------
    rho : ndarray of float, shape (lmax + 1,)
        Correlation coefficients, with bounds ``0:lmax``.
    """
    lmax = len(alm1[:,0]) - 1
    return libcurvedsky.utils.alm2rho(lmax,alm1,alm2)


def alm2cov(alm):
    r"""
    Compute auto and cross spectra between multiple alms.
  
    Parameters
    ----------
    alm : array_like of complex, shape (n, lmax + 1, lmax + 1)
        Harmonic coefficients, with bounds ``(n, 0:lmax, 0:lmax)``.
  
    Returns
    -------
    cov : ndarray of float, shape (n, n, lmax + 1)
        Auto and cross angular power spectra between ``alm[i]`` and ``alm[j]``.
    """
    
    n = len(alm[:,0,0])
    lmax = len(alm[0,:,0]) - 1

    return libcurvedsky.utils.alm2cov(n,lmax,alm)

    
def apodize(rmask,ascale,order=1,holeminsize=0):
    r"""
    Compute an apodized window function.
    
    This function partially uses Healpix's ``process_mask`` code.
    
    Parameters
    ----------
    rmask : array_like of float, shape (npix,)
        Input window function, with bounds ``0:npix-1``. Pixels with
        ``rmask == 0`` are treated as masked pixels.
    ascale : float
        Apodization length in degrees from the closest masked pixel.
    order : int, optional
        Pixel ordering. Use 1 for RING ordering; other values are treated as NESTED ordering. Default is 1.
    holeminsize : float, optional
        Minimum hole size in arcminutes. Holes smaller than this size are filled. Default is 0.
    
    Returns
    -------
    amask : ndarray of float, shape (npix,)
        Apodization window, with bounds ``0:npix-1`` and the same ordering as the input map.
    """
    
    return lib_utils.apodize(rmask,ascale,order,holeminsize)


def hp_alm2map(nside,alm,nthreads=0):
    r"""
    Transform alm coefficients to a map using Healpix ``(l, m)`` ordering.
  
    Parameters
    ----------
    nside : int
        Nside of the output map.
    alm : array_like of complex, shape (lmax + 1, mmax + 1)
        Harmonic coefficients to be transformed to a map, with bounds ``(0:lmax, 0:mmax)``.
    nthreads : int, optional
        Number of threads. Default is 0.
  
    Returns
    -------
    map : ndarray of float, shape (npix,)
        Transformed map, with bounds ``0:npix-1``.
    """

    return sht.alm2map(nside,alm,nthreads=nthreads)

    
def hp_alm2map_spin(nside,spin,elm,blm,nthreads=0):
  r"""
  Transform spin alm coefficients to maps using Healpix ``(l, m)`` ordering.
  
  Parameters
  ----------
  nside : int
      Nside of the output maps.
  spin : int
      Spin of the transform.
  elm : array_like of complex, shape (lmax + 1, mmax + 1)
      Spin-s E-like harmonic coefficients, with bounds ``(0:lmax, 0:mmax)``.
  blm : array_like of complex, shape (lmax + 1, mmax + 1)
      Spin-s B-like harmonic coefficients, with bounds ``(0:lmax, 0:mmax)``.
  nthreads : int, optional
      Number of threads. Default is 0.
  
  Returns
  -------
  map0 : ndarray of float, shape (npix,)
      Real part of the transformed map, or Q-like map, with bounds
      ``0:npix-1``.
  map1 : ndarray of float, shape (npix,)
      Imaginary part of the transformed map, or U-like map, with bounds
      ``0:npix-1``.
  """
  return sht.alm2map_spin(nside,spin,np.asarray([elm,blm]),nthreads=nthreads)

def hp_map2alm(lmax,mmax,map,nthreads=0):
  r"""
  Transform a map to alm coefficients using Healpix ``(l, m)`` ordering.
  
  Parameters
  ----------
  lmax : int
      Maximum multipole of the output alm.
  mmax : int
      Maximum azimuthal index of the output alm.
  map : array_like of float, shape (npix,)
      Input map, with bounds ``0:npix-1``.
  nthreads : int, optional
      Number of threads. Default is 0.
  
  Returns
  -------
  alm : ndarray of complex, shape (lmax + 1, mmax + 1)
      Harmonic coefficients obtained from the input map, with bounds
      ``(0:lmax, 0:mmax)``.
  """
  return sht.map2alm(lmax,mmax,map,nthreads=nthreads)

def hp_map2alm_spin(lmax,mmax,spin,map0,map1,nthreads=0):
  r"""
  Transform spin maps to alm coefficients using Healpix ``(l, m)`` ordering.
  
  For example, if ``map0`` is Q, ``map1`` is U, and ``spin == 2``, the
  output alm contains E and B modes.
  
  Parameters
  ----------
  lmax : int
      Maximum multipole of the output alm.
  mmax : int
      Maximum azimuthal index of the output alm.
  spin : int
      Spin of the transform.
  map0 : array_like of float, shape (npix,)
      Real part of the input map, with bounds ``0:npix-1``.
  map1 : array_like of float, shape (npix,)
      Imaginary part of the input map, with bounds ``0:npix-1``.
  nthreads : int, optional
      Number of threads. Default is 0.
  
  Returns
  -------
  alm : ndarray of complex, shape (2, lmax + 1, mmax + 1)
      Parity-even and parity-odd harmonic coefficients obtained from the
      input maps, with bounds ``(2, 0:lmax, 0:mmax)``.
  """
  return sht.map2alm_spin(lmax,mmax,spin,np.asarray([map0,map1]),nthreads=nthreads)

def map_mul_lfunc(imap,lfunc):
    r"""
    Filter a map by a multipole-dependent function.
    
    The input map is transformed to alm, multiplied by ``lfunc``, and then
    transformed back to a map.
    
    Parameters
    ----------
    imap : array_like of float, shape (12 * nside**2,)
        Input map, with bounds ``0:12*nside**2-1``.
    lfunc : array_like of float, shape (lmax + 1,)
        One-dimensional spectrum multiplied into the alm coefficients, with
        bounds ``0:lmax``.
    
    Returns
    -------
    omap : ndarray of float, shape (12 * nside**2,)
        Output map, with bounds ``0:12*nside**2-1``.
    """
    return lib_utils.map_mul_lfunc(imap,lfunc)

def mulwin(alm,win,nthreads=0):
    r"""
    Multiply a map reconstructed from alm by a window function.
    
    Parameters
    ----------
    alm : array_like of complex, shape (lmax + 1, mmax + 1)
        Harmonic coefficients used to construct the input map, with bounds
        ``(0:lmax, 0:mmax)``.
    win : array_like of float, shape (npix,)
        Window map, with bounds ``0:npix-1``.
    nthreads : int, optional
        Number of threads. Default is 0.
    
    Returns
    -------
    wlm : ndarray of complex, shape (lmax + 1, mmax + 1)
        Harmonic coefficients of the window-multiplied map, with bounds
        ``(0:lmax, 0:mmax)``.
    """
    
    nside = int(np.sqrt(len(win) / 12))
    lmax  = len(alm[:,0]) - 1
    mmax  = len(alm[0,:]) - 1

    map_in = sht.alm2map(nside, alm, nthreads=nthreads)
    map_w  = win * map_in
    
    return sht.map2alm(lmax, lmax, map_w, nthreads=nthreads)


def mulwin_spin(alm,win,spin=2,nthreads=0):
    r"""
    Multiply spin maps reconstructed from alm by a window function.
    
    Parameters
    ----------
    alm : array_like of complex, shape (2, lmax + 1, mmax + 1)
        Spin-s E/B-like harmonic coefficients used to construct the input
        maps, with bounds ``(2, 0:lmax, 0:mmax)``.
    win : array_like of float, shape (npix,)
        Window map, with bounds ``0:npix-1``.
    spin : int, optional
        Spin of the transform. Default is 2.
    nthreads : int, optional
        Number of threads. Default is 0.
    
    Returns
    -------
    wlm : ndarray of complex, shape (2, lmax + 1, mmax + 1)
        Parity-even and parity-odd harmonic coefficients obtained from the
        window-multiplied maps, with bounds ``(2, 0:lmax, 0:mmax)``.
    """
  
    nside = int(np.sqrt(len(win) / 12))
    lmax  = len(alm[0,:,0]) - 1
    mmax  = len(alm[0,0,:]) - 1

    map_in = sht.alm2map_spin(nside, spin, alm, nthreads=nthreads)
    map_w  = win * map_in
    
    return sht.map2alm_spin(lmax, lmax, spin, map_w, nthreads=nthreads)


def lm_healpy2healpix(almpy):
    r"""
    Convert healpy alm storage to Healpix ``(l, m)`` storage.
    
    Parameters
    ----------
    almpy : array_like of complex, shape (lmpy,)
        Healpy alm array, with bounds ``0:lmpy-1``.
    
    Returns
    -------
    almpix : ndarray of complex, shape (lmax + 1, lmax + 1)
        Healpix alm array, with bounds ``(0:lmax, 0:lmax)``.
    """
    lmpy = len(almpy)
    lmax = int((-3.+np.sqrt(1.+8.*lmpy))/2.)
    return libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax)


def lm_healpix2healpy(almpix):
    r"""
    Convert Healpix ``(l, m)`` alm storage to healpy storage.
    
    Parameters
    ----------
    almpix : array_like of complex, shape (lmax + 1, lmax + 1)
        Healpix alm array, with bounds ``(0:lmax, 0:lmax)``.
    
    Returns
    -------
    almpy : ndarray of complex, shape (lmpy,)
        Healpy alm array, with bounds ``0:lmpy-1``.
    """
    lmpy = int((lmax+1)*(lmax+2)/2.)
    lmax = len(alm[:,0]) - 1
    return libcurvedsky.utils.lm_healpy2healpix(lmax,almpix,lmpy)
    

def cosin_healpix(nside):
    r"""
    Return :math:`\cos(\theta)` as a function of the Healpix pixel index.
    
    Parameters
    ----------
    nside : int
        Nside of the desired map.
    
    Returns
    -------
    cosin : ndarray of float, shape (npix,)
        Values of :math:`\cos(\theta)`, with bounds ``0:npix-1``.
    """
    return lib_utils.cosin_healpix(nside)


def calc_mfs(nu,walm,nside=0):
    r"""
    Compute two-dimensional Minkowski functionals from a given alm.
    
    Parameters
    ----------
    nu : array_like of float, shape (bn,)
        Threshold bins, with bounds ``0:bn-1``.
    walm : array_like of complex, shape (lmax + 1, lmax + 1)
        Alm with filtering, possibly divided by the map variance, with bounds
        ``(0:lmax, 0:lmax)``.
    nside : int, optional
        Nside of the intermediate map. If 0, the closest power of 2 to
        ``3*lmax`` is used. Default is 0.
    
    Returns
    -------
    V : ndarray of float, shape (bn, 3)
        The three Minkowski functionals, ``V0``, ``V1``, and ``V2``, at each
        threshold bin, with bounds ``(0:bn-1, 0:2)``.
    """
    nu = np.asarray(nu, dtype=np.float64)

    # Determine lmax and convert alm if needed
    walm = np.asarray(walm)
    lmax = walm.shape[0] - 1
    if nside == 0:
        nside_t = 2 ** int(np.log(float(3 * lmax)) / np.log(2.0))
    else:
        nside_t = int(nside)

    der0 = sht.alm2map(nside_t, walm)
    der1 = sht.alm2map_der1(nside_t, walm)
    der2 = sht.alm2map_der2(nside_t, walm)

    return lib_utils.calc_mfs_from_derivatives(nu, der0, der1, der2)


def load_defangle_takahashi(fname,npix,verbose=False):
    r"""
    Read source-plane theta and phi coordinates from Takahashi et al. (2017).
    
    Parameters
    ----------
    fname : str
        File name.
    npix : int
        Number of pixels in the theta and phi data.
    verbose : bool, optional
        If True, print progress messages. Default is False.
    
    Returns
    -------
    theta : ndarray of float, shape (npix,)
        Theta coordinates, with bounds ``0:npix-1``.
    phi : ndarray of float, shape (npix,)
        Phi coordinates, with bounds ``0:npix-1``.
    """
    return libcurvedsky.utils.load_defangle_takahashi(fname,npix,verbose)


def polcoord2angle(theta,phi,verbose=False):
    r"""
    Convert source-plane theta and phi coordinates to deflection angles.
    
    Parameters
    ----------
    theta : array_like of float, shape (npix,)
        Theta coordinates, with bounds ``0:npix-1``.
    phi : array_like of float, shape (npix,)
        Phi coordinates, with bounds ``0:npix-1``.
    verbose : bool, optional
        If True, print progress messages. Default is False.
    
    Returns
    -------
    angle : ndarray of float, shape (2, npix)
        Deflection-angle vector containing two components, with bounds
        ``(0:1, 0:npix-1)``.
    """
    return lib_utils.polcoord2angle(theta,phi,verbose=verbose)


def polcoord2angle_alm(lmax,theta,phi,verbose=False):
    r"""
    Convert source-plane theta and phi coordinates to deflection-angle alms.
    
    Parameters
    ----------
    lmax : int
        Maximum multipole of the gradient and curl alms.
    theta : array_like of float, shape (npix,)
        Theta coordinates, with bounds ``0:npix-1``.
    phi : array_like of float, shape (npix,)
        Phi coordinates, with bounds ``0:npix-1``.
    verbose : bool, optional
        If True, print progress messages. Default is False.
    
    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Gradient-mode alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    return lib_utils.polcoord2angle_alm(lmax,theta,phi,verbose=verbose)


def mock_galaxy_takahashi(fname,zn,ngz,zi,b0=1.0,btype='sqrtz',a=2.0,b=1.0,zm=1.0,sz=0.0,zbias=0.0):
    r"""
    Compute galaxy overdensity maps from Takahashi et al. (2017) matter maps.
    
    The galaxy redshift distribution is assumed to follow Eq. (7) of
    arXiv:1810.03346.
    
    Parameters
    ----------
    fname : str
        File name of the density map.
    zn : int
        Number of tomographic bins.
    ngz : array_like of float, shape (zn,)
        Total number of galaxies in each redshift bin.
    zi : array_like of float, shape (zn + 1,)
        Redshift bin edges for the tomographic bins.
    b0 : float, optional
        Constant galaxy bias at ``z = 0``. Default is 1.0.
    btype : str, optional
        Redshift dependence type for the galaxy bias. Default is ``'sqrtz'``.
    a : float, optional
        Galaxy distribution shape parameter. Default is 2.0.
    b : float, optional
        Galaxy distribution shape parameter. Default is 1.0.
    zm : float, optional
        Mean redshift of the galaxy distribution. Default is 1.0.
    sz : float, optional
        Redshift scatter parameter. Default is 0.0.
    zbias : float, optional
        Redshift bias parameter. Default is 0.0.
    
    Returns
    -------
    gmap : ndarray of float, shape (npix, zbin)
        Galaxy number-density map for each redshift bin.
    """
    return libcurvedsky.utils.mock_galaxy_takahashi(fname,zn,ngz,zi,b0,btype,a,b,zm,sz,zbias)

