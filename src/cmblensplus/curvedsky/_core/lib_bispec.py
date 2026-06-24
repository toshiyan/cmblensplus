
import numpy as np
from .sht_ducc0 import alm2map, map2alm

PI = np.pi


#////////////////////////////////////////////////////////////////////////////////#
# Internal functions
#////////////////////////////////////////////////////////////////////////////////#

def _default_bin_pair(bp: np.ndarray):
    return np.asarray(bp[:2], dtype=np.int64)


def _middle_bin_pair(bp: np.ndarray):
    """Match Fortran bp(bn/2:bn/2+1) with 1-based indexing."""
    bn = len(bp) - 1
    start = max(bn // 2 - 1, 0)
    return np.asarray(bp[start : start + 2], dtype=np.int64)


def _resolve_sL(bp, sL=None):
    if sL is None:
        return _default_bin_pair(bp)
    arr = np.asarray(sL, dtype=np.int64)
    if arr.shape != (2,):
        raise ValueError("sL must be None or a length-2 array-like object")
    if np.sum(arr) == 0:
        return _default_bin_pair(bp)
    return arr


def _validate_bstype(bstype: str) -> str:
    bstype = str(bstype).strip().lower()
    if bstype not in {"equi", "fold", "sque", "isos"}:
        raise ValueError("bstype must be one of 'equi', 'fold', 'sque', 'isos'")
    return bstype


def _nside_from_lmax(lmax: int, bst: int) -> int:
    """Fortran: bst*2**int(log(lmax)/log(2))."""
    if lmax < 1:
        raise ValueError("lmax must be >= 1")
    return int(bst) * 2 ** int(np.log(float(lmax)) / np.log(2.0))


def _filtered_alm(alm, lmin, lmax, out_lmax = None,):
    """Return an alm array containing only lmin <= l <= lmax."""
    if out_lmax is None:
        out_lmax = lmax

    out = np.zeros((out_lmax + 1, out_lmax + 1), dtype=np.complex128)

    for ell in range(max(0, int(lmin)), min(int(lmax), out_lmax) + 1):
        out[ell, : ell + 1] = alm[ell, : ell + 1]

    return out


def _map_product_sum(kmap: np.ndarray, n: int) -> float:
    """Match the n=1,2,3 cases in the cross-bispectrum Fortran routines."""
    if n == 1:
        return float(np.sum(kmap[0] ** 3))
    if n == 2:
        return float(np.sum(kmap[0] * kmap[1] ** 2))
    if n == 3:
        return float(np.sum(kmap[0] * kmap[1] * kmap[2]))
    raise ValueError("n should satisfy 1 <= n <= 3")


#////////////////////////////////////////////////////////////////////////////////#
# Core function for bispectrum calculation
#////////////////////////////////////////////////////////////////////////////////#

def equi(lmin: int, lmax: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Equilateral shape, b[l,l,l], for a single alm."""
    nside = _nside_from_lmax(lmax, bst)
    kmap = alm2map(nside, _filtered_alm(alm, lmin, lmax), nthreads=nthreads)
    return float(np.sum(kmap ** 3) * (4.0 * PI) / kmap.size)


def fold(lmin: int, lmax: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Folded shape, b[l,l/2,l/2], for a single alm."""
    nside = _nside_from_lmax(lmax, bst)
    kmap1 = alm2map(nside, _filtered_alm(alm, lmin, lmax), nthreads=nthreads)
    kmap2 = alm2map(nside, _filtered_alm(alm, max(2, int(lmin / 2.0)), int(lmax / 2.0), out_lmax=lmax), nthreads=nthreads)
    return float(np.sum(kmap1 * kmap2 ** 2) * (4.0 * PI) / kmap1.size)


def sque(eL: np.ndarray, sL: np.ndarray, l1: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Squeezed shape, b[sL,eL,eL], for a single alm."""
    e0, e1 = map(int, eL)
    s0, s1 = map(int, sL)
    if max(s1, e1) > l1:
        raise ValueError("l1 is too small for the requested squeezed bins")
    nside = _nside_from_lmax(l1, bst)
    kmap1 = alm2map(nside, _filtered_alm(alm, s0, s1, out_lmax=l1), nthreads=nthreads)
    kmap2 = alm2map(nside, _filtered_alm(alm, e0, e1, out_lmax=l1), nthreads=nthreads)
    return float(np.sum(kmap1 * kmap2 ** 2) * (4.0 * PI) / kmap1.size)


def isos(eL: np.ndarray, aL: np.ndarray, l1: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Isosceles shape, b[eL,aL,aL], for a single alm."""
    e0, e1 = map(int, eL)
    a0, a1 = map(int, aL)
    if max(a1, e1) > l1:
        raise ValueError("l1 is too small for the requested isosceles bins")
    nside = _nside_from_lmax(l1, bst)
    kmap1 = alm2map(nside, _filtered_alm(alm, e0, e1, out_lmax=l1), nthreads=nthreads)
    kmap2 = alm2map(nside, _filtered_alm(alm, a0, a1, out_lmax=l1), nthreads=nthreads)
    return float(np.sum(kmap1 * kmap2 ** 2) * (4.0 * PI) / kmap1.size)


def xequi(lmin: int, lmax: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Cross equilateral shape for n=1,2,3 input fields."""
    n = int(alm.shape[0])
    if not 1 <= n <= 3:
        raise ValueError("n should satisfy 1 <= n <= 3")
    nside = _nside_from_lmax(lmax, bst)
    kmaps = np.asarray([
        alm2map(nside, _filtered_alm(alm[i], lmin, lmax), nthreads=nthreads)
        for i in range(n)
    ])
    return _map_product_sum(kmaps, n) * (4.0 * PI) / kmaps.shape[1]


def xfold(lmin: int, lmax: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Cross folded shape for n=1,2,3 input fields."""
    n = int(alm.shape[0])
    if not 1 <= n <= 3:
        raise ValueError("n should satisfy 1 <= n <= 3")
    nn = max(n, 2)
    nside = _nside_from_lmax(lmax, bst)
    klm = np.zeros((nn, lmax + 1, lmax + 1), dtype=np.complex128)
    for ell in range(max(0, int(lmin)), int(lmax) + 1):
        klm[0, ell, : ell + 1] = alm[0, ell, : ell + 1]
    for ell in range(max(2, int(lmin / 2.0)), int(lmax / 2.0) + 1):
        if n == 1:
            klm[1, ell, : ell + 1] = alm[0, ell, : ell + 1]
        else:
            klm[1:n, ell, : ell + 1] = alm[1:n, ell, : ell + 1]
    kmaps = np.asarray([alm2map(nside, klm[i], nthreads=nthreads) for i in range(nn)])
    return _map_product_sum(kmaps, n) * (4.0 * PI) / kmaps.shape[1]


def xsque(eL: np.ndarray, sL: np.ndarray, l1: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Cross squeezed shape for n=1,2,3 input fields."""
    n = int(alm.shape[0])
    if not 1 <= n <= 3:
        raise ValueError("n should satisfy 1 <= n <= 3")
    e0, e1 = map(int, eL)
    s0, s1 = map(int, sL)
    nn = max(n, 2)
    nside = _nside_from_lmax(l1, bst)
    klm = np.zeros((nn, l1 + 1, l1 + 1), dtype=np.complex128)
    for ell in range(s0, s1 + 1):
        klm[0, ell, : ell + 1] = alm[0, ell, : ell + 1]
    for ell in range(e0, e1 + 1):
        if n == 1:
            klm[1, ell, : ell + 1] = alm[0, ell, : ell + 1]
        else:
            klm[1:n, ell, : ell + 1] = alm[1:n, ell, : ell + 1]
    kmaps = np.asarray([alm2map(nside, klm[i], nthreads=nthreads) for i in range(nn)])
    return _map_product_sum(kmaps, n) * (4.0 * PI) / kmaps.shape[1]


def xisos(eL: np.ndarray, aL: np.ndarray, l1: int, alm: np.ndarray, bst: int = 2, nthreads: int = 0) -> float:
    """Cross isosceles shape for n=1,2,3 input fields."""
    n = int(alm.shape[0])
    if not 1 <= n <= 3:
        raise ValueError("n should satisfy 1 <= n <= 3")
    e0, e1 = map(int, eL)
    a0, a1 = map(int, aL)
    if max(a1, e1) > l1:
        raise ValueError("l1 is too small for the requested isosceles bins")
    nn = max(n, 2)
    nside = _nside_from_lmax(l1, bst)
    klm = np.zeros((nn, l1 + 1, l1 + 1), dtype=np.complex128)
    for ell in range(e0, e1 + 1):
        klm[0, ell, : ell + 1] = alm[0, ell, : ell + 1]
    for ell in range(a0, a1 + 1):
        if n == 1:
            klm[1, ell, : ell + 1] = alm[0, ell, : ell + 1]
        else:
            klm[1:n, ell, : ell + 1] = alm[1:n, ell, : ell + 1]
    kmaps = np.asarray([alm2map(nside, klm[i], nthreads=nthreads) for i in range(nn)])
    return _map_product_sum(kmaps, n) * (4.0 * PI) / kmaps.shape[1]


#////////////////////////////////////////////////////////////////////////////////#
# Interface of core bispectrum function
#////////////////////////////////////////////////////////////////////////////////#

def make_quad_gauss(alm, bst=1, nthreads=0):
    """
    Return alm for delta_NL(n) = delta_L(n) + delta_L(n)^2.

    Parameters
    ----------
    alm
        Input scalar harmonic coefficients with shape ``(lmax+1, lmax+1)``.
    nthreads
        Thread count forwarded to ``sht_ducc0``.
    """
    alm = np.asarray(alm, dtype=np.complex128)
    lmax = alm.shape[0] - 1
    if alm.shape[1] < lmax + 1:
        raise ValueError("alm must have shape (lmax+1, lmax+1)")
    nside = _nside_from_lmax(bst,lmax)
    mp = alm2map(nside, alm[: lmax + 1, : lmax + 1], nthreads=nthreads)
    mp = mp + mp ** 2
    return np.asarray(map2alm(lmax, lmax, mp, nthreads=nthreads), dtype=np.complex128)


def bispec_norm(bp, bstype="equi", bst=2, sL=None, nthreads=0):
    """Normalization of the binned reduced bispectrum for each bin."""
    bp = np.asarray(bp, dtype=np.float64)
    bn = len(bp) - 1
    if bn < 1:
        raise ValueError("bp must contain at least two bin edges")
    bstype = _validate_bstype(bstype)
    sL0 = _resolve_sL(bp, sL)
    aL = _middle_bin_pair(bp)
    norm = np.empty(bn, dtype=np.float64)

    for ib in range(bn):
        eL = np.asarray(bp[ib : ib + 2], dtype=np.int64)
        if bstype in {"equi", "fold"}:
            ilmax = int(eL[1])
        elif bstype == "sque":
            ilmax = int(max(eL[1], sL0[1]))
        else:  # isos
            ilmax = int(max(eL[1], aL[1]))

        klm = np.zeros((ilmax + 1, ilmax + 1), dtype=np.complex128)
        for ell in range(1, ilmax + 1):
            klm[ell, 0] = np.sqrt(2.0 * ell + 1.0)

        if bstype == "equi":
            norm[ib] = equi(int(eL[0]), int(eL[1]), klm, bst=bst, nthreads=nthreads)
        elif bstype == "fold":
            norm[ib] = fold(int(eL[0]), int(eL[1]), klm, bst=bst, nthreads=nthreads)
        elif bstype == "sque":
            norm[ib] = sque(eL, sL0, ilmax, klm, bst=bst, nthreads=nthreads)
        else:
            norm[ib] = isos(eL, aL, ilmax, klm, bst=bst, nthreads=nthreads)
    return norm


def bispec_bin(bp, alm, bstype="equi", bst=2, sL=None, nthreads=0):
    """Unnormalized binned reduced bispectrum for each bin."""
    bp = np.asarray(bp, dtype=np.float64)
    alm = np.asarray(alm, dtype=np.complex128)
    lmax = alm.shape[0] - 1
    bn = len(bp) - 1
    if bp[-1] > lmax:
        raise ValueError("not enough size of alm")
    bstype = _validate_bstype(bstype)
    sL0 = _resolve_sL(bp, sL)
    aL = _middle_bin_pair(bp)
    bis = np.empty(bn, dtype=np.float64)

    for ib in range(bn):
        eL = np.asarray(bp[ib : ib + 2], dtype=np.int64)
        if bstype in {"equi", "fold"}:
            ilmax = int(eL[1])
        elif bstype == "sque":
            ilmax = int(max(eL[1], sL0[1]))
        else:
            ilmax = int(max(eL[1], aL[1]))

        alm_tmp = alm[: ilmax + 1, : ilmax + 1]
        if bstype == "equi":
            bis[ib] = equi(int(eL[0]), int(eL[1]), alm_tmp, bst=bst, nthreads=nthreads)
        elif bstype == "fold":
            bis[ib] = fold(int(eL[0]), int(eL[1]), alm_tmp, bst=bst, nthreads=nthreads)
        elif bstype == "sque":
            bis[ib] = sque(eL, sL0, ilmax, alm[: ilmax + 1, : ilmax + 1], bst=bst, nthreads=nthreads)
        else:
            bis[ib] = isos(eL, aL, ilmax, alm[: ilmax + 1, : ilmax + 1], bst=bst, nthreads=nthreads)
    return bis


def xbispec_bin(bp, alm, bstype="equi", bst=2, sL=None, nthreads=0):
    """Unnormalized binned reduced cross-bispectrum for n=1,2,3 fields."""
    bp = np.asarray(bp, dtype=np.float64)
    alm = np.asarray(alm, dtype=np.complex128)
    n = alm.shape[0]
    if not 1 <= n <= 3:
        raise ValueError("n should satisfy 1 <= n <= 3")
    lmax = alm.shape[1] - 1
    bn = len(bp) - 1
    if bp[-1] > lmax:
        raise ValueError("not enough size of alm")
    bstype = _validate_bstype(bstype)
    sL0 = _resolve_sL(bp, sL)
    aL = _middle_bin_pair(bp)
    bis = np.empty(bn, dtype=np.float64)

    for ib in range(bn):
        eL = np.asarray(bp[ib : ib + 2], dtype=np.int64)
        if bstype == "equi":
            bis[ib] = xequi(int(eL[0]), int(eL[1]), alm[:, : int(eL[1]) + 1, : int(eL[1]) + 1], bst=bst, nthreads=nthreads)
        elif bstype == "fold":
            bis[ib] = xfold(int(eL[0]), int(eL[1]), alm[:, : int(eL[1]) + 1, : int(eL[1]) + 1], bst=bst, nthreads=nthreads)
        elif bstype == "sque":
            l1 = int(max(eL[1], sL0[1]))
            bis[ib] = xsque(eL, sL0, l1, alm[:, : l1 + 1, : l1 + 1], bst=bst, nthreads=nthreads)
        else:
            l1 = int(max(eL[1], aL[1]))
            bis[ib] = xisos(eL, aL, l1, alm[:, : l1 + 1, : l1 + 1], bst=bst, nthreads=nthreads)
    return bis


    
