import numpy as np
from .sht_ducc0 import alm2map_spin, map2alm


def _default_nside(lmax: int) -> int:
    if int(lmax) < 1:
        return 1
    return 2 ** int(np.log(float(lmax)) / np.log(2.0))


def _validate_square_alm(name: str, alm: np.ndarray, lmax: int) -> np.ndarray:
    arr = np.asarray(alm, dtype=np.complex128)
    if arr.ndim != 2 or arr.shape[0] < lmax + 1 or arr.shape[1] < lmax + 1:
        raise ValueError(f"{name} must have shape at least ({lmax + 1}, {lmax + 1})")
    return arr


def _validate_spectrum(name: str, spec: np.ndarray, lmax: int) -> np.ndarray:
    arr = np.asarray(spec, dtype=np.float64)
    if arr.ndim != 1 or arr.shape[0] < lmax + 1:
        raise ValueError(f"{name} must have length at least {lmax + 1}")
    return arr


def qeb(
    lmax: int,
    rlmin: int,
    rlmax: int,
    EB: np.ndarray,
    Elm: np.ndarray,
    Blm: np.ndarray,
    nside_t: int = 0,
    verbose: bool = False,
    nthreads: int = 0,
):
    """
    Odd EB quadratic estimator for amplitude modulation.
    """
    lmax = int(lmax)
    rlmin = int(rlmin)
    rlmax = int(rlmax)
    nside = int(nside_t) if int(nside_t) != 0 else _default_nside(lmax)
    if verbose:
        print(f"calc tau-EB odd estimator with nside= {nside}")

    EB = _validate_spectrum("EB", EB, rlmax)
    Elm = _validate_square_alm("Elm", Elm, rlmax)
    Blm = _validate_square_alm("Blm", Blm, rlmax)

    # convolution 1
    xlm = np.zeros((2, rlmax + 1, rlmax + 1), dtype=np.complex128)
    for ell in range(max(0, rlmin), rlmax + 1):
        xlm[0, ell, : ell + 1] = Elm[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, xlm, nthreads=nthreads), dtype=np.float64)

    xlm.fill(0.0)
    for ell in range(max(0, rlmin), rlmax + 1):
        xlm[1, ell, : ell + 1] = EB[ell] * Blm[ell, : ell + 1]
    A2 = np.asarray(alm2map_spin(nside, 2, xlm, nthreads=nthreads), dtype=np.float64)

    mp = -A[0] * A2[1] + A[1] * A2[0]

    # convolution 2
    xlm.fill(0.0)
    for ell in range(max(0, rlmin), rlmax + 1):
        xlm[0, ell, : ell + 1] = EB[ell] * Elm[ell, : ell + 1]
    A2 = np.asarray(alm2map_spin(nside, 2, xlm, nthreads=nthreads), dtype=np.float64)

    xlm.fill(0.0)
    for ell in range(max(0, rlmin), rlmax + 1):
        xlm[1, ell, : ell + 1] = Blm[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, xlm, nthreads=nthreads), dtype=np.float64)

    mp = mp + A[0] * A2[1] - A[1] * A2[0]

    xlm_scalar = np.asarray(map2alm(lmax, lmax, mp, nthreads=nthreads), dtype=np.complex128)

    alm = np.zeros((lmax + 1, lmax + 1), dtype=np.complex128)
    for ell in range(1, lmax + 1):
        alm[ell, : ell + 1] = -2.0 * xlm_scalar[ell, : ell + 1]
    return alm
