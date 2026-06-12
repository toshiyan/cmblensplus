
import numpy as np
from dataclasses import dataclass, field
from . import sht_ducc0 as sht
from .. import libcurvedsky # wrapped fortran sources


@dataclass
class MgCovmat:
    nij: np.ndarray | None = None       # [n, pix]
    clh: np.ndarray | None = None       # [n, n, l]
    nl: np.ndarray | None = None        # [n, l, m]
    imap: np.ndarray | None = None      # [n, pix]

@dataclass
class MgChain:
    ytype: str
    n: int
    mnmax: np.ndarray
    nside: np.ndarray
    lmax: np.ndarray
    itn: np.ndarray
    eps: np.ndarray
    verbose: bool = False
    ro: int = 50
    cv: list[list[MgCovmat]] = field(default_factory=list)
    lsp: int = 0
    matn: int = 0
    pmat: np.ndarray | None = None
    mat: np.ndarray | None = None
    imat: np.ndarray | None = None
    npix: np.ndarray | None = None


def set_mgchain(ytype: str, chn: int, mn: int, mnmaxs, lmaxs, nsides, itns, eps, verbose=False, ro=50) -> MgChain:
    mnmaxs = np.asarray(mnmaxs, dtype=int)
    lmaxs = np.asarray(lmaxs, dtype=int)
    nsides = np.asarray(nsides, dtype=int)
    itns = np.asarray(itns, dtype=int)
    eps = np.asarray(eps, dtype=float)
    if nsides.shape != (chn, mn):
        raise ValueError(f"nsides must have shape ({chn}, {mn}), got {nsides.shape}")
    mgc = MgChain(ytype=ytype, n=chn, mnmax=mnmaxs, nside=nsides, lmax=lmaxs, itn=itns, eps=eps, verbose=verbose, ro=ro)
    mgc.npix = 12 * nsides**2
    mgc.lsp = int(lmaxs[-1])
    mgc.cv = [[MgCovmat() for _ in range(mn)] for _ in range(chn)]
    if verbose:
        print("nside:", mgc.nside)
        print("#maps:", mgc.mnmax)
        print("lmax: ", mgc.lmax)
        print("#iter:", mgc.itn)
        print("eps:  ", mgc.eps)
        print(mgc.lsp)
        if mgc.ro != 50:
            print(f"residual output per {mgc.ro} iters")
    return mgc


def cg_algorithm(
    n: int, 
    lmax: int, 
    b: np.ndarray, 
    mgc: MgChain, 
    chain: int = 0, 
    ratio: np.ndarray | None = None
) -> np.ndarray:
    """
    Conjugate-gradient solver for A x = b. chain is zero-based.
    """
    lmaxch = int(mgc.lmax[chain])
    if lmax != lmaxch:
        raise ValueError("lmax is inconsistent")
    b = np.asarray(b, dtype=np.complex128)
    absb = np.sqrt(np.sum(np.abs(b) ** 2))
    if not np.isfinite(absb):
        raise FloatingPointError("|b|^2 is NaN")

    power = 1 if mgc.ytype == "scal" else 2

    if chain == 0 and mgc.pmat is None:
        mgc.pmat = np.zeros((n, lmax + 1, lmax + 1), dtype=float)
        for ni in range(n):
            mm = np.zeros((lmax + 1, lmax + 1), dtype=float)
            for mi in range(int(mgc.mnmax[0])):
                cv = mgc.cv[0][mi]
                ave_nij = float(np.mean(cv.nij[ni]))
                if cv.nl is not None:
                    ave_nl = np.mean(cv.nl[ni, :, 0])
                else:
                    ave_nl = 1.0
                for ell in range(lmax + 1):
                    mm[ell, : ell + 1] += (cv.clh[ni, ni, ell] ** power) * ave_nij * ave_nl
            mgc.pmat[ni] = 1.0 / (1.0 + mm)

    x = precondition(n, lmax, b, mgc, chain)
    if not np.all(np.isfinite(x)):
        raise FloatingPointError("Mb is NaN")

    r = b - matmul_lhs(mgc.ytype, n, int(mgc.mnmax[chain]), lmax, mgc.cv[chain][: int(mgc.mnmax[chain])], x)
    z = precondition(n, lmax, r, mgc, chain)
    p = z.copy()
    d0 = np.vdot(r, z)
    d = d0

    for i in range(int(mgc.itn[chain])):
        ap = matmul_lhs(mgc.ytype, n, int(mgc.mnmax[chain]), lmax, mgc.cv[chain][: int(mgc.mnmax[chain])], p)
        denom = np.vdot(p, ap)
        alpha = d / denom
        x = x + alpha * p
        r = r - alpha * ap
        absr = np.sqrt(np.sum(np.abs(r) ** 2))
        if chain == 0 and ratio is not None:
            ratio[i] = absr / absb
        if chain == 0 and mgc.verbose and (i + 1) % mgc.ro == 0:
            print(i + 1, absr / absb, d / d0)
        z = precondition(n, lmax, r, mgc, chain)
        td = np.vdot(r, z)
        p = z + (td / d) * p
        d = td
        if absr < mgc.eps[chain] * absb:
            if chain == 0 and mgc.verbose:
                print(i + 1, absr / absb)
            break
    return x


def precondition(n: int, lmax: int, r: np.ndarray, mgc: MgChain, chain: int) -> np.ndarray:
    """
    Precondition matrix for conjugate-gradient decent
    """
    if mgc.n == 1:
        return mgc.pmat[:, : lmax + 1, : lmax + 1] * r
    x = np.zeros_like(r)
    lmax0 = int(mgc.lmax[chain + 1])
    if chain + 1 == mgc.n - 1:
        x[:, : lmax0 + 1, : lmax0 + 1] = libcurvedsky.utils.densemat_multi(r[:, : lmax0 + 1, : lmax0 + 1], mgc.imat, lmax0)
    else:
        x[:, : lmax0 + 1, : lmax0 + 1] = cg_algorithm(n, lmax0, r[:, : lmax0 + 1, : lmax0 + 1], mgc, chain + 1)
    if lmax0 + 1 <= int(mgc.lmax[chain]):
        x[:, lmax0 + 1 : lmax + 1, :] = mgc.pmat[:, lmax0 + 1 : lmax + 1, : lmax + 1] * r[:, lmax0 + 1 : lmax + 1, :]
    return x


def matmul_rhs(ytype: str, n: int, mn: int, lmax: int, cv_list: list[MgCovmat]) -> np.ndarray:
    rhs = np.zeros((n, lmax + 1, lmax + 1), dtype=np.complex128)
    for mi in range(mn):
        if ytype == "scal":
            rhs += matmul_rhs_scal(n, lmax, cv_list[mi])
        elif ytype == "cmb":
            rhs += matmul_rhs_cmb(n, lmax, cv_list[mi])
        else:
            raise ValueError("no correct ytype assigned")
    return rhs


def _nside_from_map(map_in: np.ndarray) -> int:
    return int(np.sqrt(map_in.shape[-1] / 12.0))


def _map2alm_components(map_in: np.ndarray, lmax: int) -> np.ndarray:
    """Use sht.map2alm. n=3 is treated as T plus spin-2 Q/U."""
    map_in = np.asarray(map_in, dtype=float)
    n = map_in.shape[0]
    if n == 3:
        t = sht.map2alm(lmax, lmax, map_in[0])
        eb = sht.map2alm_spin(lmax, lmax, 2, map_in[1:3])
        return np.asarray([t, eb[0], eb[1]], dtype=np.complex128)
    if n == 2:
        return np.asarray(sht.map2alm_spin(lmax, lmax, 2, map_in), dtype=np.complex128)
    return np.asarray([sht.map2alm(lmax, lmax, map_in[0])], dtype=np.complex128)


def _alm2map_components(alm: np.ndarray, nside: int) -> np.ndarray:
    """Use sht.alm2map. n=3 is treated as T plus spin-2 E/B."""
    alm = np.asarray(alm, dtype=np.complex128)
    n = alm.shape[0]
    if n == 3:
        t = sht.alm2map(nside, alm[0])
        qu = sht.alm2map_spin(nside, 2, alm[1:3])
        return np.asarray([t, qu[0], qu[1]], dtype=float)
    if n == 2:
        return np.asarray(sht.alm2map_spin(nside, 2, alm), dtype=float)
    return np.asarray([sht.alm2map(nside, alm[0])], dtype=float)


def matmul_rhs_cmb(n: int, lmax: int, cv: MgCovmat) -> np.ndarray:
    mp = np.asarray(cv.imap, dtype=float) * np.asarray(cv.nij, dtype=float)
    if cv.nl is not None:
        blm = _map2alm_components(mp, lmax)
        blm = blm * cv.nl
        mp = _alm2map_components(blm, _nside_from_map(mp))
        mp = mp * cv.nij
    blm = _map2alm_components(mp, lmax)
    return matmul_cov_alm(cv.clh, blm)


def matmul_rhs_scal(n: int, lmax: int, cv: MgCovmat) -> np.ndarray:
    alm = np.zeros((n, lmax + 1, lmax + 1), dtype=np.complex128)
    nside = _nside_from_map(cv.nij)
    for ni in range(n):
        mp = np.asarray([cv.imap[ni] * cv.nij[ni]], dtype=float)
        if cv.nl is not None:
            blm = np.asarray([sht.map2alm(lmax, lmax, mp[0])], dtype=np.complex128)
            blm[0] *= cv.nl[ni]
            mp = np.asarray([sht.alm2map(nside, blm[0])], dtype=float)
            mp[0] *= cv.nij[ni]
        alm[ni] = sht.map2alm(lmax, lmax, mp[0])
    return alm


def matmul_lhs(ytype: str, n: int, mn: int, lmax: int, cv_list: list[MgCovmat], x: np.ndarray) -> np.ndarray:
    v = np.asarray(x, dtype=np.complex128).copy()
    for mi in range(mn):
        if ytype == "scal":
            v += matmul_lhs_scal(n, lmax, cv_list[mi], x)
        elif ytype == "cmb":
            v += matmul_lhs_cmb(n, lmax, cv_list[mi], x)
        else:
            raise ValueError("no correct ytype assigned")
    return v


def matmul_lhs_cmb(n: int, lmax: int, cv: MgCovmat, x: np.ndarray) -> np.ndarray:
    alm = matmul_cov_alm(cv.clh, x)
    mp = _alm2map_components(alm, _nside_from_map(cv.nij))
    mp = mp * cv.nij
    if cv.nl is not None:
        alm = _map2alm_components(mp, lmax)
        alm = alm * cv.nl
        mp = _alm2map_components(alm, _nside_from_map(cv.nij))
        mp = mp * cv.nij
    alm = _map2alm_components(mp, lmax)
    return matmul_cov_alm(cv.clh, alm)


def matmul_lhs_scal(n: int, lmax: int, cv: MgCovmat, x: np.ndarray) -> np.ndarray:
    cov_alm = matmul_cov_alm(cv.clh, x)
    v = np.zeros((n, lmax + 1, lmax + 1), dtype=np.complex128)
    nside = _nside_from_map(cv.nij)
    for ni in range(n):
        blm = cov_alm[ni]
        mp = np.asarray([sht.alm2map(nside, blm)], dtype=float)
        mp[0] *= cv.nij[ni]
        if cv.nl is not None:
            blm = sht.map2alm(lmax, lmax, mp[0])
            blm = blm * cv.nl[ni]
            mp = np.asarray([sht.alm2map(nside, blm)], dtype=float)
            mp[0] *= cv.nij[ni]
        v[ni] = sht.map2alm(lmax, lmax, mp[0])
    return v


def matmul_cov_alm(cov: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Compute cov[n,n,l] @ x[n,l,m] for each (l,m)."""
    cov = np.asarray(cov, dtype=float)
    x = np.asarray(x, dtype=np.complex128)
    n, _, lmax1 = cov.shape
    out = np.zeros((n, lmax1, lmax1), dtype=np.complex128)
    for ell in range(lmax1):
        out[:, ell, : ell + 1] = cov[:, :, ell] @ x[:, ell, : ell + 1]
    return out


def coarse_invmatrix(mgc: MgChain, n: int, lmax: int) -> None:
    mn = int(mgc.mnmax[-1])
    mgc.matn = n * (lmax + 1) * (lmax + 2) // 2
    a0 = np.zeros((mgc.matn, mgc.matn), dtype=np.complex128)
    i = 0
    for ni in range(n):
        for ell in range(lmax + 1):
            for m in range(ell + 1):
                x = np.zeros((n, lmax + 1, lmax + 1), dtype=np.complex128)
                x[ni, ell, m] = 1.0
                mx = matmul_lhs(mgc.ytype, n, mn, lmax, mgc.cv[-1][:mn], x)
                col0 = libcurvedsky.utils.trans_alm2array_1d(mx, mgc.matn)
                if m != 0:
                    x.fill(0.0)
                    x[ni, ell, m] = 1j
                    mx = matmul_lhs(mgc.ytype, n, mn, lmax, mgc.cv[-1][:mn], x)
                    col1 = libcurvedsky.utils.trans_alm2array_1d(mx, mgc.matn)
                    col0 = (col0 + col1 / 1j) / 2.0
                a0[:, i] = col0
                i += 1
    mgc.imat = np.linalg.inv(a0)


