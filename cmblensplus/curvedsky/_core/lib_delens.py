
import numpy as np
from .sht_ducc0 import alm2map_spin, map2alm_spin
try:
    from .sht_ducc0 import alm2map_der
except ImportError:
    alm2map_der = None


def _default_nside(lmax: int) -> int:
    """Fortran default: ``2**int(log(lmax)/log(2))``."""
    if lmax < 1:
        return 1
    return 2 ** int(np.log(float(lmax)) / np.log(2.0))


def _validate_square_alm(name: str, alm: np.ndarray, lmax: int) -> np.ndarray:
    arr = np.asarray(alm, dtype=np.complex128)
    if arr.ndim != 2 or arr.shape[0] < lmax + 1 or arr.shape[1] < lmax + 1:
        raise ValueError(f"{name} must have shape at least ({lmax + 1}, {lmax + 1})")
    return arr


def lensingb(
    lmax: int,
    elmin: int,
    elmax: int,
    plmin: int,
    plmax: int,
    wElm: np.ndarray,
    wplm: np.ndarray,
    nside_t: int = 0,
    gtype: str = "p",
    nthreads: int = 0,
) -> np.ndarray:
    """
    Compute lensing B-mode alm from filtered E-mode and lensing potential.
    """
    lmax = int(lmax)
    elmin = int(elmin)
    elmax = int(elmax)
    plmin = int(plmin)
    plmax = int(plmax)

    wElm = _validate_square_alm("wElm", wElm, elmax)
    wplm = _validate_square_alm("wplm", wplm, plmax)

    nside = int(nside_t)
    if nside == 0:
        nside = _default_nside(max(elmax, plmax))

    gtype = str(gtype).strip().lower()[:1] or "p"
    if gtype not in {"p", "k"}:
        raise ValueError("gtype must be 'p' for phi or 'k' for kappa")

    ilk = np.ones(plmax + 1, dtype=np.float64)
    if gtype == "k":
        for ell in range(1, plmax + 1):
            ilk[ell] = 2.0 / float(ell * (ell + 1))

    # A1: spin-1 transform of sqrt((l+2)(l-1)/2) * E_lm
    alm_e = np.zeros((2, elmax + 1, elmax + 1), dtype=np.complex128)
    for ell in range(max(elmin, 0), elmax + 1):
        fac = np.sqrt(float((ell + 2) * (ell - 1)) * 0.5)
        alm_e[0, ell, : ell + 1] = wElm[ell, : ell + 1] * fac
    A1 = np.asarray(alm2map_spin(nside, 1, alm_e, nthreads=nthreads), dtype=np.float64)

    # A3: spin-3 transform of sqrt((l-2)(l+3)/2) * E_lm
    alm_e.fill(0.0)
    for ell in range(max(elmin, 0), elmax + 1):
        fac_arg = float((ell - 2) * (ell + 3)) * 0.5
        fac = np.sqrt(fac_arg) if fac_arg >= 0.0 else 0.0
        alm_e[0, ell, : ell + 1] = wElm[ell, : ell + 1] * fac
    A3 = np.asarray(alm2map_spin(nside, 3, alm_e, nthreads=nthreads), dtype=np.float64)

    # A: spin-1 transform of sqrt(l(l+1)/2) * phi_lm.
    # If the input is kappa, multiply by 2/[l(l+1)]
    alm_p = np.zeros((2, plmax + 1, plmax + 1), dtype=np.complex128)
    for ell in range(max(plmin, 0), plmax + 1):
        fac = np.sqrt(float(ell * (ell + 1)) * 0.5) * ilk[ell]
        alm_p[0, ell, : ell + 1] = wplm[ell, : ell + 1] * fac
    A = np.asarray(alm2map_spin(nside, 1, alm_p, nthreads=nthreads), dtype=np.float64)

    if A.shape[0] != 2 or A1.shape[0] != 2 or A3.shape[0] != 2:
        raise ValueError("alm2map_spin must return arrays with shape (2, npix)")

    # convolution
    conv = np.empty_like(A)
    conv[0] = A[0] * (A1[0] - A3[0]) - A[1] * (A1[1] + A3[1])
    conv[1] = A[0] * (A1[1] - A3[1]) + A[1] * (A1[0] + A3[0])

    alm_b = np.asarray(map2alm_spin(lmax, lmax, 2, conv, nthreads=nthreads), dtype=np.complex128)
    return alm_b[1, : lmax + 1, : lmax + 1]



def phi2grad(nside: int, lmax: int, plm: np.ndarray) -> FloatArray:
    """Return the deflection vector at each Healpix pixel.

    This is the Python equivalent of the f2py-facing Fortran routine where
    the ``npix`` argument is replaced by ``nside``.

    Returns an array with shape ``(12*nside**2, 2)``.
    """
    nside = int(nside)
    lmax = int(lmax)
    plm = _validate_square_alm("plm", plm, lmax)
    _mp, alpha, _dalpha = _call_alm2map_der(nside, lmax, plm)
    return np.asarray(alpha, dtype=np.float64)


def shiftvec(nside: int, lmax: int, plm: np.ndarray, nremap: int = 3) -> FloatArray:
    """Return the anti-deflection vector ``beta`` for iterative delensing.

    The vector satisfies approximately ``beta(n) + alpha(n + beta(n)) = 0``.
    This mirrors the Fortran iteration order exactly.

    Returns an array with shape ``(12*nside**2, 2)``.
    """
    nside = int(nside)
    lmax = int(lmax)
    nremap = int(nremap)
    plm = _validate_square_alm("plm", plm, lmax)

    _mp, alpha, dalpha = _call_alm2map_der(nside, lmax, plm)
    if dalpha is None:
        raise NotImplementedError("shiftvec() requires alm2map_der to return second derivatives as dalpha")

    alpha = np.asarray(alpha, dtype=np.float64)
    dalpha = np.asarray(dalpha, dtype=np.float64)
    npix = 12 * nside ** 2
    if alpha.shape != (npix, 2):
        raise ValueError("alm2map_der must return alpha with shape (npix, 2)")
    if dalpha.shape != (npix, 3):
        raise ValueError("alm2map_der must return dalpha with shape (npix, 3)")

    beta = np.zeros((npix, 2), dtype=np.float64)
    for _ in range(nremap):
        beta[:, 0] = alpha[:, 0] - dalpha[:, 0] * beta[:, 0] - dalpha[:, 1] * beta[:, 1]
        beta[:, 1] = alpha[:, 1] - dalpha[:, 1] * beta[:, 0] - dalpha[:, 2] * beta[:, 1]

    return -beta
