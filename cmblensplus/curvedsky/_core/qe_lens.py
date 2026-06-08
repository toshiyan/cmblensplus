from . import utils
import numpy


def _validate_alm(name: str, alm: Array, rlmax: int) -> Array:
    arr = np.asarray(alm, dtype=np.complex128)
    if arr.shape[0] < rlmax + 1 or arr.shape[1] < rlmax + 1:
        raise ValueError(f"{name} must have shape at least ({rlmax + 1}, {rlmax + 1}); got {arr.shape}")
    return arr


def _validate_cl(name: str, cl: Array, rlmax: int) -> Array:
    arr = np.asarray(cl, dtype=np.float64)
    if arr.shape[0] < rlmax + 1:
        raise ValueError(f"{name} must have length at least {rlmax + 1}; got {arr.shape[0]}")
    return arr
    
def _default_nside(lmax: int, nside_t: int) -> int:
    if nside_t != 0:
        return int(nside_t)
    if lmax <= 0:
        return 1
    # Fortran used: 2**int(log(lmax)/log(2)).
    return 2 ** int(np.log(float(lmax)) / np.log(2.0))


def _ilk(lmax: int, gtype: str) -> Array:
    out = np.ones(lmax + 1, dtype=np.float64)
    if gtype == "k":
        ell = np.arange(1, lmax + 1, dtype=np.float64)
        out[1:] = 2.0 / (ell * (ell + 1.0))
    return out


def _zeros_spin(lmax: int, mmax: Optional[int] = None) -> Array:
    if mmax is None:
        mmax = lmax
    return np.zeros((2, lmax + 1, mmax + 1), dtype=np.complex128)


def _zeros_scalar(lmax: int, mmax: Optional[int] = None) -> Array:
    if mmax is None:
        mmax = lmax
    return np.zeros((1, lmax + 1, mmax + 1), dtype=np.complex128)


def _fill_component(alm: Array, component: int, source: Array, rlmin: int, rlmax: int, weights: Optional[Array] = None) -> None:
    for ell in range(rlmin, rlmax + 1):
        if weights is None:
            alm[component, ell, : ell + 1] = source[ell, : ell + 1]
        else:
            alm[component, ell, : ell + 1] = weights[ell] * source[ell, : ell + 1]


def _finish_gradient_curl(lmax: int, spin_alm: Array, ilk: Array, prefactor: float = 1.0) -> Tuple[Array, Array]:
    glm = np.zeros((lmax + 1, lmax + 1), dtype=np.complex128)
    clm = np.zeros((lmax + 1, lmax + 1), dtype=np.complex128)
    for ell in range(1, lmax + 1):
        fac = ilk[ell] * prefactor * np.sqrt(float(ell * (ell + 1)))
        glm[ell, : ell + 1] = fac * spin_alm[0, ell, : ell + 1]
        clm[ell, : ell + 1] = fac * spin_alm[1, ell, : ell + 1]
    return glm, clm


def _common_spin_product(A: Array, A1: Array, A3: Array) -> Array:
    """
    map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
    map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
    """

    maps = np.empty_like(A, dtype=np.float64)
    maps[:, 0] = A[:, 0] * (A1[:, 0] - A3[:, 0]) + A[:, 1] * (A1[:, 1] - A3[:, 1])
    maps[:, 1] = -A[:, 0] * (A1[:, 1] + A3[:, 1]) + A[:, 1] * (A1[:, 0] + A3[:, 0])
    return maps


def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t=0,gtype='',verbose=False):
    """
    Reconstructing CMB lensing potential and its curl mode from the temperature quadratic estimator

    Args:
        :lmax (int): Maximum multipole of output lensing potential alms
        :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
        :fC [l] (double): TT spectrum, with bounds (0:rlmax)
        :Tlm1 [l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
        :Tlm2 [l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

    Args(optional):
        :nside_t (int): Nside for the convolution calculation
        :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
        :verbose (bool): Output messages, default to False

    Returns:
        :glm [l,m] (dcmplx): CMB lensing potential alm, with bounds (0:lmax,0:lmax)
        :clm [l,m] (dcmplx): Curl mode (pseudo lensing potential) alm, with bounds (0:lmax,0:lmax)

    """

    backend = _require_backend()
    fC = _validate_cl("fC", fC, rlmax)
    Tlm1 = _validate_alm("Tlm1", Tlm1, rlmax)
    Tlm2 = _validate_alm("Tlm2", Tlm2, rlmax)

    nside = _default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qTT lens estimator with nside= {nside}")

    ilk = _ilk(lmax, gtype)

    # convolution
    alm1 = _zeros_scalar(rlmax)
    _fill_component(alm1, 0, Tlm1, rlmin, rlmax)
    at = utils.alm2map(nside, rlmax, rlmax, alm1)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    weights = fC * np.sqrt(ell * (ell + 1.0))
    alm1 = _zeros_spin(rlmax)
    _fill_component(alm1, 0, Tlm2, rlmin, rlmax, weights)
    maps = utils.alm2map_spin(nside, rlmax, rlmax, 1, alm1)
    maps[:, 0] *= at
    maps[:, 1] *= at
    blm = utils.map2alm_spin(nside, lmax, lmax, 1, maps)

    # compute glm and clm
    return _finish_gradient_curl(lmax, blm, ilk, prefactor=1.0)


def qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the TE quadratic estimator

  Args:
    :lmax (int): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
    :fC [l] (double): TE spectrum, with bounds (0:rlmax)
    :Tlm [l,m] (dcmplx): Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Elm [l,m] (dcmplx): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (int): Nside for the convolution calculation
    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (bool): Output messages, default to False

  Returns:
    :glm [l,m] (dcmplx): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [l,m] (dcmplx): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  """
  return libcurvedsky.rec_lens.qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside_t,gtype,verbose)

def qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    :lmax (int): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
    :fC [l] (double): TE spectrum, with bounds (0:rlmax)
    :Tlm [l,m] (dcmplx): Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Blm [l,m] (dcmplx): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (int): Nside for the convolution calculation
    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (bool): Output messages, default to False

  Returns:
    :glm [l,m] (dcmplx): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [l,m] (dcmplx): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  """
  return libcurvedsky.rec_lens.qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside_t,gtype,verbose)

def qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the EE quadratic estimator

  Args:
    :lmax (int): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
    :fC [l] (double): EE spectrum, with bounds (0:rlmax)
    :Elm1 [l,m] (dcmplx): 1st inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Elm2 [l,m] (dcmplx): 2nd inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (int): Nside for the convolution calculation
    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (bool): Output messages, default to False

  Returns:
    :glm [l,m] (dcmplx): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [l,m] (dcmplx): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  """
  return libcurvedsky.rec_lens.qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside_t,gtype,verbose)

def qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :lmax (int): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
    :fC [l] (double): EE spectrum, with bounds (0:rlmax)
    :Elm [l,m] (dcmplx): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [l,m] (dcmplx): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (int): Nside for the convolution calculation
    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (bool): Output messages, default to False

  Returns:
    :glm [l,m] (dcmplx): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [l,m] (dcmplx): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  """
  return libcurvedsky.rec_lens.qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside_t,gtype,verbose)

def qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the BB quadratic estimator

  Args:
    :lmax (int): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (int): Minimum/Maximum multipoles of CMB for reconstruction
    :fC [l] (double): BB spectrum, with bounds (0:rlmax)
    :Blm1 [l,m] (dcmplx): 1st inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm2 [l,m] (dcmplx): 2nd inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (int): Nside for the convolution calculation
    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (bool): Output messages, default to False

  Returns:
    :glm [l,m] (dcmplx): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [l,m] (dcmplx): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  """
  return libcurvedsky.rec_lens.qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside_t,gtype,verbose)

