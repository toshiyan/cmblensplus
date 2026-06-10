import numpy as np
from .sht_ducc0 import alm2map, alm2map_spin, map2alm_spin


def _default_nside(lmax: int) -> int:
    """Fortran default: ``2**int(log(lmax)/log(2))``."""
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


def _ilk(lmax: int, gtype: str) -> np.ndarray:
    ilk = np.ones(lmax + 1, dtype=np.float64)
    if (str(gtype).strip().lower()[:1] == "k"):
        for ell in range(1, lmax + 1):
            ilk[ell] = 2.0 / float(ell * (ell + 1))
    return ilk


def _finish_lensing_outputs(
    lmax: int,
    gtype: str,
    zlm: np.ndarray,
) -> tuple[ComplexArray, ComplexArray]:
    ilk = _ilk(lmax, gtype)
    glm = np.zeros((lmax + 1, lmax + 1), dtype=np.complex128)
    clm = np.zeros((lmax + 1, lmax + 1), dtype=np.complex128)
    for ell in range(1, lmax + 1):
        fac = ilk[ell] * np.sqrt(float(ell * (ell + 1)))
        glm[ell, : ell + 1] = fac * zlm[0, ell, : ell + 1]
        clm[ell, : ell + 1] = fac * zlm[1, ell, : ell + 1]
    return glm, clm


def qte(
    lmax: int,
    rlmin: int,
    rlmax: int,
    fC: np.ndarray,
    Tlm: np.ndarray,
    Elm: np.ndarray,
    nside_t: int = 0,
    gtype: str = "",
    verbose: bool = False,
    nthreads: int = 0,
):
    """Imaginary lensing TE quadratic estimator. Returns ``(glm, clm)``."""
    lmax = int(lmax)
    rlmin = int(rlmin)
    rlmax = int(rlmax)
    nside = int(nside_t) if int(nside_t) != 0 else _default_nside(lmax)
    if verbose:
        print(f"calc qTE ilens estimator with nside= {nside}")

    fC = _validate_spectrum("fC", fC, rlmax)
    Tlm = _validate_square_alm("Tlm", Tlm, rlmax)
    Elm = _validate_square_alm("Elm", Elm, rlmax)

    alm_e = np.zeros((2, rlmax + 1, rlmax + 1), dtype=np.complex128)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm_e[0, ell, : ell + 1] = Elm[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, alm_e, nthreads=nthreads), dtype=np.float64)

    alm1 = np.zeros_like(alm_e)
    alm3 = np.zeros_like(alm_e)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm1[0, ell, : ell + 1] = fC[ell] * Tlm[ell, : ell + 1] * np.sqrt(float((ell + 2) * (ell - 1)))
        alm3[0, ell, : ell + 1] = fC[ell] * Tlm[ell, : ell + 1] * np.sqrt(float((ell - 2) * (ell + 3)))
    A1 = np.asarray(alm2map_spin(nside, 1, alm1, nthreads=nthreads), dtype=np.float64)
    A3 = np.asarray(alm2map_spin(nside, 3, alm3, nthreads=nthreads), dtype=np.float64)

    mp = np.empty_like(A)
    mp[0] = -A[0] * (A1[1] - A3[1]) + A[1] * (A1[0] - A3[0])
    mp[1] = -A[0] * (A1[0] + A3[0]) - A[1] * (A1[1] + A3[1])

    zlm = np.asarray(map2alm_spin(lmax, lmax, 1, mp, nthreads=nthreads), dtype=np.complex128)
    return _finish_lensing_outputs(lmax, gtype, zlm)


def qtb(
    lmax: int,
    rlmin: int,
    rlmax: int,
    fC: np.ndarray,
    Tlm: np.ndarray,
    Blm: np.ndarray,
    nside_t: int = 0,
    gtype: str = "",
    verbose: bool = False,
    nthreads: int = 0,
):
    """Imaginary lensing TB quadratic estimator. Returns ``(glm, clm)``."""
    lmax = int(lmax)
    rlmin = int(rlmin)
    rlmax = int(rlmax)
    nside = int(nside_t) if int(nside_t) != 0 else _default_nside(lmax)
    if verbose:
        print(f"calc qTB ilens estimator with nside= {nside}")

    fC = _validate_spectrum("fC", fC, rlmax)
    Tlm = _validate_square_alm("Tlm", Tlm, rlmax)
    Blm = _validate_square_alm("Blm", Blm, rlmax)

    alm = np.zeros((2, rlmax + 1, rlmax + 1), dtype=np.complex128)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm[1, ell, : ell + 1] = Blm[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, alm, nthreads=nthreads), dtype=np.float64)

    alm1 = np.zeros_like(alm)
    alm3 = np.zeros_like(alm)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm1[0, ell, : ell + 1] = fC[ell] * Tlm[ell, : ell + 1] * np.sqrt(float((ell + 2) * (ell - 1)))
        alm3[0, ell, : ell + 1] = fC[ell] * Tlm[ell, : ell + 1] * np.sqrt(float((ell - 2) * (ell + 3)))
    A1 = np.asarray(alm2map_spin(nside, 1, alm1, nthreads=nthreads), dtype=np.float64)
    A3 = np.asarray(alm2map_spin(nside, 3, alm3, nthreads=nthreads), dtype=np.float64)

    tlm = np.zeros((rlmax + 1, rlmax + 1), dtype=np.complex128)
    for ell in range(max(0, rlmin), rlmax + 1):
        tlm[ell, : ell + 1] = Tlm[ell, : ell + 1]
    AT = np.asarray(alm2map(nside, tlm, nthreads=nthreads), dtype=np.float64)

    alm_ab = np.zeros_like(alm)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm_ab[1, ell, : ell + 1] = Blm[ell, : ell + 1] * fC[ell] * np.sqrt(float(ell * (ell + 1)))
    AB = np.asarray(alm2map_spin(nside, 1, alm_ab, nthreads=nthreads), dtype=np.float64)

    mp = np.empty_like(A)
    mp[0] = -A[0] * (A1[1] - A3[1]) + A[1] * (A1[0] - A3[0]) + 2.0 * AT * AB[1]
    mp[1] = -A[0] * (A1[0] + A3[0]) - A[1] * (A1[1] + A3[1]) - 2.0 * AT * AB[0]

    zlm = np.asarray(map2alm_spin(lmax, lmax, 1, mp, nthreads=nthreads), dtype=np.complex128)
    return _finish_lensing_outputs(lmax, gtype, zlm)


def qee(
    lmax: int,
    rlmin: int,
    rlmax: int,
    fC: np.ndarray,
    Elm1: np.ndarray,
    Elm2: np.ndarray,
    nside_t: int = 0,
    gtype: str = "",
    verbose: bool = False,
    nthreads: int = 0,
):
    """Imaginary lensing EE quadratic estimator. Returns ``(glm, clm)``."""
    lmax = int(lmax)
    rlmin = int(rlmin)
    rlmax = int(rlmax)
    nside = int(nside_t) if int(nside_t) != 0 else _default_nside(lmax)
    if verbose:
        print(f"calc qEE ilens estimator with nside= {nside}")

    fC = _validate_spectrum("fC", fC, rlmax)
    Elm1 = _validate_square_alm("Elm1", Elm1, rlmax)
    Elm2 = _validate_square_alm("Elm2", Elm2, rlmax)

    alm = np.zeros((2, rlmax + 1, rlmax + 1), dtype=np.complex128)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm[0, ell, : ell + 1] = Elm1[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, alm, nthreads=nthreads), dtype=np.float64)

    alm1 = np.zeros_like(alm)
    blm = np.zeros_like(alm)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm1[0, ell, : ell + 1] = fC[ell] * Elm2[ell, : ell + 1] * np.sqrt(float((ell + 2) * (ell - 1)))
        # Preserves the factor used in rec_ilens.f90: sqrt((l-1)*(l+3)).
        blm[0, ell, : ell + 1] = fC[ell] * Elm2[ell, : ell + 1] * np.sqrt(float((ell - 1) * (ell + 3)))
    A1 = np.asarray(alm2map_spin(nside, 1, alm1, nthreads=nthreads), dtype=np.float64)
    A3 = np.asarray(alm2map_spin(nside, 3, blm, nthreads=nthreads), dtype=np.float64)

    mp = np.empty_like(A)
    mp[0] = -A[0] * (A1[1] - A3[1]) + A[1] * (A1[0] - A3[0])
    mp[1] = -A[0] * (A1[0] + A3[0]) - A[1] * (A1[1] + A3[1])

    zlm = np.asarray(map2alm_spin(lmax, lmax, 1, mp, nthreads=nthreads), dtype=np.complex128)
    return _finish_lensing_outputs(lmax, gtype, zlm)


def qeb(
    lmax: int,
    rlmin: int,
    rlmax: int,
    fC: np.ndarray,
    Elm: np.ndarray,
    Blm: np.ndarray,
    nside_t: int = 0,
    gtype: str = "",
    verbose: bool = False,
    nthreads: int = 0,
):
    """Imaginary lensing EB quadratic estimator. Returns ``(glm, clm)``."""
    lmax = int(lmax)
    rlmin = int(rlmin)
    rlmax = int(rlmax)
    nside = int(nside_t) if int(nside_t) != 0 else _default_nside(lmax)
    if verbose:
        print(f"calc qEB ilens estimator with nside= {nside}")

    fC = _validate_spectrum("fC", fC, rlmax)
    Elm = _validate_square_alm("Elm", Elm, rlmax)
    Blm = _validate_square_alm("Blm", Blm, rlmax)

    alm = np.zeros((2, rlmax + 1, rlmax + 1), dtype=np.complex128)

    # first part
    for ell in range(max(0, rlmin), rlmax + 1):
        alm[1, ell, : ell + 1] = Blm[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, alm, nthreads=nthreads), dtype=np.float64)

    alm1 = np.zeros_like(alm)
    alm3 = np.zeros_like(alm)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm1[0, ell, : ell + 1] = fC[ell] * Elm[ell, : ell + 1] * np.sqrt(float((ell + 2) * (ell - 1)))
        alm3[0, ell, : ell + 1] = fC[ell] * Elm[ell, : ell + 1] * np.sqrt(float((ell - 2) * (ell + 3)))
    A1 = np.asarray(alm2map_spin(nside, 1, alm1, nthreads=nthreads), dtype=np.float64)
    A3 = np.asarray(alm2map_spin(nside, 3, alm3, nthreads=nthreads), dtype=np.float64)

    mp = np.empty_like(A)
    mp[0] = -A[0] * (A1[1] - A3[1]) + A[1] * (A1[0] - A3[0])
    mp[1] = -A[0] * (A1[0] + A3[0]) - A[1] * (A1[1] + A3[1])

    # second part
    alm.fill(0.0)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm[0, ell, : ell + 1] = Elm[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, alm, nthreads=nthreads), dtype=np.float64)

    alm1.fill(0.0)
    alm3.fill(0.0)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm1[1, ell, : ell + 1] = fC[ell] * Blm[ell, : ell + 1] * np.sqrt(float((ell + 2) * (ell - 1)))
        alm3[1, ell, : ell + 1] = fC[ell] * Blm[ell, : ell + 1] * np.sqrt(float((ell - 2) * (ell + 3)))
    A1 = np.asarray(alm2map_spin(nside, 1, alm1, nthreads=nthreads), dtype=np.float64)
    A3 = np.asarray(alm2map_spin(nside, 3, alm3, nthreads=nthreads), dtype=np.float64)

    mp[0] = mp[0] + A[0] * (A1[1] - A3[1]) - A[1] * (A1[0] - A3[0])
    mp[1] = mp[1] + A[0] * (A1[0] + A3[0]) + A[1] * (A1[1] + A3[1])

    zlm = np.asarray(map2alm_spin(lmax, lmax, 1, mp, nthreads=nthreads), dtype=np.complex128)
    return _finish_lensing_outputs(lmax, gtype, zlm)


def qbb(
    lmax: int,
    rlmin: int,
    rlmax: int,
    fC: np.ndarray,
    Blm1: np.ndarray,
    Blm2: np.ndarray,
    nside_t: int = 0,
    gtype: str = "",
    verbose: bool = False,
    nthreads: int = 0,
):
    """Imaginary lensing BB quadratic estimator. Returns ``(glm, clm)``."""
    lmax = int(lmax)
    rlmin = int(rlmin)
    rlmax = int(rlmax)
    nside = int(nside_t) if int(nside_t) != 0 else _default_nside(lmax)
    if verbose:
        print(f"calc qBB ilens estimator with nside= {nside}")

    fC = _validate_spectrum("fC", fC, rlmax)
    Blm1 = _validate_square_alm("Blm1", Blm1, rlmax)
    Blm2 = _validate_square_alm("Blm2", Blm2, rlmax)

    alm = np.zeros((2, rlmax + 1, rlmax + 1), dtype=np.complex128)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm[1, ell, : ell + 1] = Blm1[ell, : ell + 1]
    A = np.asarray(alm2map_spin(nside, 2, alm, nthreads=nthreads), dtype=np.float64)

    alm1 = np.zeros_like(alm)
    alm3 = np.zeros_like(alm)
    for ell in range(max(0, rlmin), rlmax + 1):
        alm1[1, ell, : ell + 1] = fC[ell] * Blm2[ell, : ell + 1] * np.sqrt(float((ell + 2) * (ell - 1)))
        alm3[1, ell, : ell + 1] = fC[ell] * Blm2[ell, : ell + 1] * np.sqrt(float((ell - 2) * (ell + 3)))
    A1 = np.asarray(alm2map_spin(nside, 1, alm1, nthreads=nthreads), dtype=np.float64)
    A3 = np.asarray(alm2map_spin(nside, 3, alm3, nthreads=nthreads), dtype=np.float64)

    mp = np.empty_like(A)
    mp[0] = A[0] * (A1[1] - A3[1]) - A[1] * (A1[0] - A3[0])
    mp[1] = A[0] * (A1[0] + A3[0]) + A[1] * (A1[1] + A3[1])

    zlm = np.asarray(map2alm_spin(lmax, lmax, 1, mp, nthreads=nthreads), dtype=np.complex128)
    return _finish_lensing_outputs(lmax, gtype, zlm)
