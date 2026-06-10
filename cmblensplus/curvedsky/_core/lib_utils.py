
import numpy as np


#//// utilities for quadratic estimators ////#

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


def _zeros_scal(lmax: int, mmax: Optional[int] = None) -> Array:
    if mmax is None:
        mmax = lmax
    return np.zeros((1, lmax + 1, mmax + 1), dtype=np.complex128)


def _fill_component(alm_type, component: int, source: Array, rlmin: int, rlmax: int, weights: Optional[Array] = None) -> None:
    """
    component = 0 if source is Talm or Ealm
    component = 1 if source is Balm
    """

    if alm_type == 'scal':
        alm = _zeros_scal(rlmax)
    if alm_type == 'spin':
        alm = _zeros_spin(rlmax)

    for ell in range(rlmin, rlmax + 1):
        if weights is None:
            alm[component, ell, : ell + 1] = source[ell, : ell + 1]
        else:
            alm[component, ell, : ell + 1] = weights[ell] * source[ell, : ell + 1]
    return alm


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
    maps[0] =  A[0] * (A1[0] - A3[0]) + A[1] * (A1[1] - A3[1])
    maps[1] = -A[0] * (A1[1] + A3[1]) + A[1] * (A1[0] + A3[0])
    return maps



#//// utilities for C^-1 ////#


def lmpy2lmax(lmpy: int) -> int:
    """
    Return lmax from the length of a healpy-packed alm array.
    """
    return int((-3.0 + np.sqrt(1.0 + 8.0 * int(lmpy))) / 2.0)


def lmax_to_nside(lmax: int) -> int:
    """
    Return the nearest upper power-of-two nside for roughly lmax / 3.
    """
    nside_real = float(lmax) / 3.0
    power_of_two = 0
    while 2**power_of_two < nside_real:
        power_of_two += 1
    return 2**power_of_two


def alm_size(n: int, lmax: int) -> int:
    return int(n * (lmax + 1) * (lmax + 2) // 2)


def lm_indices(lmax: int):
    """
    Return indices (ell, m) for valid alm entries with 0 <= m <= ell <= lmax.
    """
    return np.tril_indices(lmax + 1)

    
def trans_alm2array_1d(alm: np.ndarray, arrn: int | None = None) -> np.ndarray:
    """
    Pack alm[n,l,m] with 0 <= m <= l into a 1D array.
    """
    alm = np.asarray(alm)
    n = alm.shape[0]
    lmax = alm.shape[1] - 1

    if arrn is None:
        arrn = alm_size(n, lmax)

    required = alm_size(n, lmax)
    if arrn < required:
        raise ValueError(f"wrong size: arrn={arrn}, required={required}")

    ell, m = lm_indices(lmax)

    arr = np.zeros(arrn, dtype=np.complex128)
    arr[:required] = alm[:, ell, m].reshape(-1)

    return arr


def trans_array2alm_1d(arr: np.ndarray, n: int, lmax: int) -> np.ndarray:
    """
    Unpack a 1D array into alm[n,l,m] with entries only for m <= l.
    """
    arr = np.asarray(arr, dtype=np.complex128)

    required = alm_size(n, lmax)
    if len(arr) < required:
        raise ValueError(f"wrong size: len(arr)={len(arr)}, required={required}")

    ell, m = lm_indices(lmax)

    alm = np.zeros((n, lmax + 1, lmax + 1), dtype=np.complex128)
    alm[:, ell, m] = arr[:required].reshape(n, -1)

    return alm


def trans_array2alm_2d(arr: np.ndarray, n: int, lmax: int) -> np.ndarray:
    """
    Unpack a dense matrix into alm[n,l,m,n,l,m].
    """
    arr = np.asarray(arr, dtype=np.complex128)

    required = alm_size(n, lmax)
    if arr.shape[0] < required or arr.shape[1] < required:
        raise ValueError(
            f"wrong size: arr.shape={arr.shape}, required={required}x{required}"
        )

    ell, m = lm_indices(lmax)
    nlm = len(ell)

    alm = np.zeros(
        (n, lmax + 1, lmax + 1, n, lmax + 1, lmax + 1),
        dtype=np.complex128,
    )

    block = arr[:required, :required].reshape(n, nlm, n, nlm)

    alm[:, ell, m, :, ell, m] = block

    return alm


