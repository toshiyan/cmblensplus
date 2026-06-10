
from __future__ import annotations

import time
import numpy as np
import healpy as hp
from .hp_cgd import ( MgCovmat, set_mgchain, clhalf, coarse_invmatrix, cg_algorithm, correct_filtering, matmul_rhs, )


def _nside_from_npix(npix: int) -> int:
    nside = int(round(np.sqrt(int(npix) / 12.0)))
    if 12 * nside**2 != int(npix):
        raise ValueError(f"npix={npix} is not a valid HEALPix pixel count")
    return nside


def ud_grade_invnoise(map_in: np.ndarray, nside_out: int, order_in: str = "RING", order_out: str = "RING") -> np.ndarray:
    """Downgrade/upgrade an inverse-noise or hit-count-like map.

    ``power=-2`` makes downgrade conserve the sum over subpixels, matching the
    usual treatment of inverse-variance weights / hit maps rather than simple
    averaging.
    """
    return np.asarray(
        hp.ud_grade(np.asarray(map_in, dtype=float), nside_out=nside_out,
                    order_in=order_in, order_out=order_out, power=-2),
        dtype=float,
    )


def _as_1d_param(x, chn: int, default, dtype=float, name: str = "param") -> np.ndarray:
    if x is None:
        arr = np.full(chn, default, dtype=dtype)
    else:
        arr = np.asarray(x, dtype=dtype)
        if arr.size == 1 and chn > 1:
            arr = np.full(chn, arr.item(), dtype=dtype)
    if arr.shape != (chn,):
        raise ValueError(f"{name} must have length {chn}, got shape {arr.shape}")
    return arr


def _fill_nl(mgc, inl: np.ndarray, by_frequency: bool) -> None:
    """Populate cv.nl when an inverse-noise spectrum is provided."""
    if np.sum(inl) == 0:
        return
    for c in range(mgc.n):
        lmaxc = int(mgc.lmax[c])
        for mi in range(len(mgc.cv[c])):
            cv = mgc.cv[c][mi]
            n = cv.nij.shape[0]
            nl = np.zeros((n, lmaxc + 1, lmaxc + 1), dtype=float)
            for ni in range(n):
                for ell in range(2, lmaxc + 1):
                    if by_frequency:
                        nl[ni, ell, : ell + 1] = inl[ni, mi, ell]
                    else:
                        nl[ni, ell, : ell + 1] = inl[ni, ell]
            cv.nl = nl


def cnfilter_freq(cl: np.ndarray, bl: np.ndarray, iNcov: np.ndarray, maps: np.ndarray,
    *, chn: int = 1, lmaxs=None, nsides=None, itns=None, eps=None, filter: str = "W", inl=None, verbose: bool = False, ro: int = 50, stat: str = "",
                 ) -> np.ndarray:
    """
    Combine multiple frequency CMB maps optimally.
    """
    cl = np.asarray(cl, dtype=float)
    bl = np.asarray(bl, dtype=float)
    iNcov = np.asarray(iNcov, dtype=float)
    maps = np.asarray(maps, dtype=float)
    n = cl.shape[0]
    mn = bl.shape[0]
    lmax = cl.shape[1] - 1
    npix = maps.shape[2]
    if inl is None:
        inl = np.zeros((n, mn, lmax + 1), dtype=float)
    else:
        inl = np.asarray(inl, dtype=float)

    eps = _as_1d_param(eps, chn, 1e-6, float, "eps")
    itns = _as_1d_param(itns, chn, 1, int, "itns")
    if chn == 1:
        ilmaxs = np.array([lmax], dtype=int)
        insides = np.full((1, mn), _nside_from_npix(npix), dtype=int)
        mnmaxs = np.array([mn], dtype=int)
    else:
        lmaxs = _as_1d_param(lmaxs, chn, 0, int, "lmaxs")
        nsides = _as_1d_param(nsides, chn, 0, int, "nsides")
        if lmax != int(lmaxs[0]):
            raise ValueError(f"input lmax is wrong: {lmax},{lmaxs[0]}")
        if npix != 12 * int(nsides[0]) ** 2:
            raise ValueError(f"input npix0 is wrong: {npix},{12 * int(nsides[0]) ** 2}")
        ilmaxs = lmaxs
        insides = np.repeat(nsides[:, None], mn, axis=1)
        mnmaxs = np.full(chn, mn, dtype=int)
        mnmaxs[1:] = 1

    clh = clhalf(cl, bl)
    mgc = set_mgchain("cmb", chn, mn, mnmaxs, ilmaxs, insides, itns, eps, verbose, ro)

    for mi in range(mn):
        mgc.cv[0][mi].imap = maps[:, mi, :].copy()
        for c in range(mgc.n):
            mgc.cv[c][mi].nij = np.zeros((n, int(mgc.npix[c, mi])), dtype=float)
            mgc.cv[c][mi].clh = np.zeros((n, n, int(mgc.lmax[c]) + 1), dtype=float)
        if iNcov.shape[2] != int(mgc.npix[0, mi]):
            raise ValueError(f"iNcov size is wrong: {iNcov.shape[2]},{mgc.npix[0, mi]}")
        mgc.cv[0][mi].clh = clh[:, :, mi, :].copy()
        mgc.cv[0][mi].nij = iNcov[:, mi, :].copy()

    _fill_nl(mgc, inl, by_frequency=True)

    b = matmul_rhs(mgc.ytype, n, mn, lmax, mgc.cv[0])
    for mi in range(mn):
        mgc.cv[0][mi].imap = None

    if mgc.n > 1:
        if stat or verbose:
            print("degrade inv noise cov")
        for c in range(1, mgc.n):
            for mi in range(mn):
                for ni in range(n):
                    mgc.cv[c][mi].nij[ni] = ud_grade_invnoise(mgc.cv[0][mi].nij[ni], int(mgc.nside[c, mi]))
                if mi == 0:
                    continue
                if mgc.cv[c][mi].nij.shape[1] != int(mgc.npix[c, 0]):
                    raise ValueError("iNcov size is inconsistent")
                mgc.cv[c][0].nij += mgc.cv[c][mi].nij
            mgc.cv[c][0].clh = np.mean(clh[:, :, :, : int(mgc.lmax[c]) + 1], axis=2)
        coarse_invmatrix(mgc, n, int(mgc.lmax[-1]))

    t1 = time.time()
    xlm = cg_algorithm(n, lmax, b, mgc, 0)
    if stat or verbose:
        print("real time:", time.time() - t1)
    return correct_filtering(cl, xlm, filter)


def cnfilter_kappa(cov: np.ndarray, iNcov: np.ndarray, maps: np.ndarray,
    *, chn: int = 1, lmaxs=None, nsides=None, itns=None, eps=None, inl=None, verbose: bool = False, ro: int = 50, stat: str = "",
                  ) -> np.ndarray:
    """
    Compute inverse-variance weighted multipoles for kappa maps.
    """
    cov = np.asarray(cov, dtype=float)
    iNcov = np.asarray(iNcov, dtype=float)
    maps = np.asarray(maps, dtype=float)
    n = cov.shape[0]
    lmax = cov.shape[2] - 1
    npix = maps.shape[1]
    if inl is None:
        inl = np.zeros((n, lmax + 1), dtype=float)
    else:
        inl = np.asarray(inl, dtype=float)

    eps = _as_1d_param(eps, chn, 1e-6, float, "eps")
    itns = _as_1d_param(itns, chn, 1, int, "itns")
    if chn == 1:
        ilmaxs = np.array([lmax], dtype=int)
        insides = np.full((1, 1), _nside_from_npix(npix), dtype=int)
    else:
        lmaxs = _as_1d_param(lmaxs, chn, 0, int, "lmaxs")
        nsides = _as_1d_param(nsides, chn, 0, int, "nsides")
        if lmax != int(lmaxs[0]):
            raise ValueError(f"input lmax is wrong: {lmax},{lmaxs[0]}")
        if npix != 12 * int(nsides[0]) ** 2:
            raise ValueError(f"input npix0 is wrong: {npix},{12 * int(nsides[0]) ** 2}")
        ilmaxs = lmaxs
        insides = nsides[:, None]
    mnmaxs = np.ones(chn, dtype=int)
    mgc = set_mgchain("scal", chn, 1, mnmaxs, ilmaxs, insides, itns, eps, verbose, ro)

    mgc.cv[0][0].imap = maps.copy()
    for c in range(mgc.n):
        mgc.cv[c][0].nij = np.zeros((n, int(mgc.npix[c, 0])), dtype=float)
        mgc.cv[c][0].clh = np.zeros((n, n, int(mgc.lmax[c]) + 1), dtype=float)
    if iNcov.shape[1] != int(mgc.npix[0, 0]):
        raise ValueError(f"iNcov size is wrong: {iNcov.shape[1]},{mgc.npix[0, 0]}")
    mgc.cv[0][0].clh = cov.copy()
    mgc.cv[0][0].nij = iNcov.copy()

    _fill_nl(mgc, inl, by_frequency=False)

    b = matmul_rhs(mgc.ytype, n, 1, lmax, mgc.cv[0])
    mgc.cv[0][0].imap = None

    if mgc.n > 1:
        if stat or verbose:
            print("degrade inv noise cov")
        for c in range(1, mgc.n):
            for ni in range(n):
                mgc.cv[c][0].nij[ni] = ud_grade_invnoise(mgc.cv[0][0].nij[ni], int(mgc.nside[c, 0]))
            mgc.cv[c][0].clh = cov[:, :, : int(mgc.lmax[c]) + 1]
        coarse_invmatrix(mgc, n, int(mgc.lmax[-1]))

    t1 = time.time()
    xlm = cg_algorithm(n, lmax, b, mgc, 0)
    if stat or verbose:
        print("real time:", time.time() - t1)
    return xlm


def cnfilter_freq_nside(cl: np.ndarray, bl0: np.ndarray, bl1: np.ndarray, iNcov0: np.ndarray, iNcov1: np.ndarray, maps0: np.ndarray, maps1: np.ndarray,
    *, chn=1, lmaxs=None, nsides0=None, nsides1=None, itns=None, eps=None, filter="W", inl=None, verbose=False, reducmn=0, ro=50, stat="",
                       ) -> np.ndarray:
    """
    Same as cnfilter_freq, but with two input Nsides.
    """
    cl = np.asarray(cl, dtype=float)
    bl0 = np.asarray(bl0, dtype=float)
    bl1 = np.asarray(bl1, dtype=float)
    iNcov0 = np.asarray(iNcov0, dtype=float)
    iNcov1 = np.asarray(iNcov1, dtype=float)
    maps0 = np.asarray(maps0, dtype=float)
    maps1 = np.asarray(maps1, dtype=float)
    n = cl.shape[0]
    mn0 = bl0.shape[0]
    mn1 = bl1.shape[0]
    mn = mn0 + mn1
    lmax = cl.shape[1] - 1
    npix0 = maps0.shape[2]
    npix1 = maps1.shape[2]
    if inl is None:
        inl = np.zeros((n, mn, lmax + 1), dtype=float)
    else:
        inl = np.asarray(inl, dtype=float)

    bl = np.concatenate([bl0, bl1], axis=0)
    clh = clhalf(cl, bl)

    eps = _as_1d_param(eps, chn, 1e-6, float, "eps")
    itns = _as_1d_param(itns, chn, 1, int, "itns")
    mnmaxs = np.full(chn, mn, dtype=int)
    if chn == 1:
        ilmaxs = np.array([lmax], dtype=int)
        insides = np.empty((1, mn), dtype=int)
        insides[0, :mn0] = _nside_from_npix(npix0)
        insides[0, mn0:] = _nside_from_npix(npix1)
    else:
        lmaxs = _as_1d_param(lmaxs, chn, 0, int, "lmaxs")
        nsides0 = _as_1d_param(nsides0, chn, 0, int, "nsides0")
        nsides1 = _as_1d_param(nsides1, chn, 0, int, "nsides1")
        if lmax != int(lmaxs[0]):
            raise ValueError(f"input lmax is wrong: {lmax},{lmaxs[0]}")
        if npix0 != 12 * int(nsides0[0]) ** 2:
            raise ValueError(f"input npix0 is wrong: {npix0},{12 * int(nsides0[0]) ** 2}")
        if npix1 != 12 * int(nsides1[0]) ** 2:
            raise ValueError(f"input npix1 is wrong: {npix1},{12 * int(nsides1[0]) ** 2}")
        ilmaxs = lmaxs
        insides = np.empty((chn, mn), dtype=int)
        for c in range(chn):
            insides[c, :mn0] = nsides0[c]
            insides[c, mn0:] = nsides1[c]
        if reducmn == 1:
            mnmaxs[1:] = 2
            insides[1:, 1] = nsides1[1:]
        elif reducmn == 2:
            if chn >= 2:
                mnmaxs[1] = 2
                insides[1, 1] = nsides1[1]
            if chn >= 3:
                mnmaxs[2:] = 1

    mgc = set_mgchain("cmb", chn, mn, mnmaxs, ilmaxs, insides, itns, eps, verbose, ro)

    for mi in range(mn):
        mgc.cv[0][mi].imap = (maps0[:, mi, :] if mi < mn0 else maps1[:, mi - mn0, :]).copy()
        for c in range(mgc.n):
            mgc.cv[c][mi].nij = np.zeros((n, int(mgc.npix[c, mi])), dtype=float)
            mgc.cv[c][mi].clh = np.zeros((n, n, int(mgc.lmax[c]) + 1), dtype=float)
        mgc.cv[0][mi].clh = clh[:, :, mi, :].copy()
        if mi < mn0:
            if iNcov0.shape[2] != int(mgc.npix[0, mi]):
                raise ValueError(f"iNcov0 size is wrong: {iNcov0.shape[2]},{mgc.npix[0, mi]}")
            mgc.cv[0][mi].nij = iNcov0[:, mi, :].copy()
        else:
            if iNcov1.shape[2] != int(mgc.npix[0, mi]):
                raise ValueError(f"iNcov1 size is wrong: {iNcov1.shape[2]},{mgc.npix[0, mi]}")
            mgc.cv[0][mi].nij = iNcov1[:, mi - mn0, :].copy()

    _fill_nl(mgc, inl, by_frequency=True)

    b = matmul_rhs(mgc.ytype, n, mn, lmax, mgc.cv[0])
    for mi in range(mn):
        mgc.cv[0][mi].imap = None

    if mgc.n > 1:
        if stat or verbose:
            print("degrade inv noise cov")
        for c in range(1, mgc.n):
            mm = int(mgc.mnmax[c])
            if mm >= 4:
                if mm != mn:
                    raise ValueError("Inconsistency in the number of maps in the multigrid chain.")
                for mi in range(mn):
                    mgc.cv[c][mi].clh = clh[:, :, mi, : int(mgc.lmax[c]) + 1]
                    for ni in range(n):
                        mgc.cv[c][mi].nij[ni] = ud_grade_invnoise(mgc.cv[0][mi].nij[ni], int(mgc.nside[c, mi]))
            elif mm == 2:
                for rn in range(2):
                    nside = int(nsides0[c] if rn == 0 else nsides1[c])
                    mgc.cv[c][rn].nij.fill(0.0)
                    start = 0 if rn == 0 else mn0
                    stop = mn0 if rn == 0 else mn0 + mn1
                    for mi in range(start, stop):
                        nij = np.zeros((n, 12 * nside**2), dtype=float)
                        for ni in range(n):
                            nij[ni] = ud_grade_invnoise(mgc.cv[0][mi].nij[ni], nside)
                        mgc.cv[c][rn].nij += nij
                mgc.cv[c][0].clh = np.mean(clh[:, :, :mn0, : int(mgc.lmax[c]) + 1], axis=2)
                mgc.cv[c][1].clh = np.mean(clh[:, :, mn0:, : int(mgc.lmax[c]) + 1], axis=2)
            elif mm == 1:
                nside = int(mgc.nside[c, 0])
                mgc.cv[c][0].nij.fill(0.0)
                for mi in range(mn):
                    nij = np.zeros((n, 12 * nside**2), dtype=float)
                    for ni in range(n):
                        nij[ni] = ud_grade_invnoise(mgc.cv[0][mi].nij[ni], nside)
                    mgc.cv[c][0].nij += nij
                mgc.cv[c][0].clh = np.mean(clh[:, :, :, : int(mgc.lmax[c]) + 1], axis=2)
            else:
                raise ValueError("Specified multigrid chain is not supported. Please change reducmn.")
        coarse_invmatrix(mgc, n, int(mgc.lmax[-1]))

    t1 = time.time()
    xlm = cg_algorithm(n, lmax, b, mgc, 0)
    if stat or verbose:
        print("real time:", time.time() - t1)
    return correct_filtering(cl, xlm, filter)
