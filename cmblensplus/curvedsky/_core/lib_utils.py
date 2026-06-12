
import numpy as np
import healpy as hp
from scipy.spatial import cKDTree
from collections import deque


def nside_from_npix(npix: int) -> int:
    nside = int(round(np.sqrt(npix / 12.0)))
    if 12 * nside * nside != npix:
        raise ValueError(f"npix={npix} is not a valid HEALPix npix = 12*nside**2")
    return nside


def default_apod_window(s, a):
    s = np.asarray(s, dtype=np.float64)

    s_out = 1.0
    s_in = s_out * a

    w = np.zeros_like(s, dtype=np.float64)

    inside = s < s_in
    trans = (s >= s_in) & (s <= s_out)

    w[inside] = 1.0

    if np.any(trans):
        x = (s_out - s[trans]) / (s_out - s_in)
        w[trans] = x - np.sin(2.0 * np.pi * x) / (2.0 * np.pi)

    if w.ndim == 0:
        return float(w)
    return w

    
def subpatch_mask(nside_full: int, nside_sub: int, ipix_sub: int, ascale: float, apod_window=default_apod_window):
    """
    Based on the fortran code written by Ryo Nagata.
    """
    npix_full = 12*nside_full**2

    ipix = np.arange(npix_full, dtype=np.int64)

    # Full-map pixel coordinates
    theta, phi = hp.pix2ang(nside_full, ipix, nest=False)

    # Identify subpatch for every full-resolution pixel
    patch_id = hp.ang2pix(nside_sub, theta, phi, nest=False)

    mask = np.zeros(npix_full, dtype=np.float64)
    in_patch = patch_id == ipix_sub

    if not np.any(in_patch):
        return mask

    # Center of selected subpatch
    theta_cent, phi_cent = hp.pix2ang(nside_sub, ipix_sub, nest=False)

    th = theta[in_patch]
    ph = phi[in_patch]

    cos_dtheta = (
        np.cos(th) * np.cos(theta_cent)
        + np.sin(th) * np.cos(ph) * np.sin(theta_cent) * np.cos(phi_cent)
        + np.sin(th) * np.sin(ph) * np.sin(theta_cent) * np.sin(phi_cent)
    )

    # Numerical safety for arccos
    cos_dtheta = np.clip(cos_dtheta, -1.0, 1.0)
    dtheta = np.arccos(cos_dtheta)

    mask[in_patch] = apod_window(dtheta, ascale)

    return mask


def get_baseline(npix: int, nside_subpatch: int, QU):
    """
    Based on the fortran code written by Ryo Nagata.
    """
    QU = np.asarray(QU, dtype=np.float64)
    if QU.shape != (npix, 2):
        raise ValueError(f"QU must have shape ({npix}, 2), got {QU.shape}")

    nside = nside_from_npix(npix)
    npix_subpatch = 12 * nside_subpatch**2

    ipix = np.arange(npix, dtype=np.int64)
    theta, phi = hp.pix2ang(nside, ipix, nest=False)

    patch_map = hp.ang2pix(nside_subpatch, theta, phi, nest=False)

    counts = np.bincount(patch_map, minlength=npix_subpatch).astype(np.float64)

    qsum = np.bincount(patch_map, weights=QU[:, 0], minlength=npix_subpatch)
    usum = np.bincount(patch_map, weights=QU[:, 1], minlength=npix_subpatch)

    patch_pmean = np.zeros((npix_subpatch, 2), dtype=np.float64)

    nonzero = counts > 0
    patch_pmean[nonzero, 0] = qsum[nonzero] / counts[nonzero]
    patch_pmean[nonzero, 1] = usum[nonzero] / counts[nonzero]

    blmap = patch_pmean[patch_map]

    return blmap

def get_winmap(nside_large: int, nside_small: int, ipix_pix: int, apod: float, apod_window=default_apod_window) -> float:
    """
    Based on the fortran code written by Ryo Nagata.
    """
    theta_pix, phi_pix = hp.pix2ang(nside_small, ipix_pix, nest=False)

    ipix_sub = hp.ang2pix(nside_large, theta_pix, phi_pix, nest=False)
    theta_sub, phi_sub = hp.pix2ang(nside_large, ipix_sub, nest=False)

    cos_dtheta = (
        np.cos(theta_pix) * np.cos(theta_sub)
        + np.sin(theta_pix) * np.cos(phi_pix) * np.sin(theta_sub) * np.cos(phi_sub)
        + np.sin(theta_pix) * np.sin(phi_pix) * np.sin(theta_sub) * np.sin(phi_sub)
    )

    cos_dtheta = np.clip(cos_dtheta, -1.0, 1.0)
    dtheta = np.arccos(cos_dtheta)

    return float(apod_window(np.array([dtheta]), apod)[0])





def eb_separate(lmax, W, Q, U):
    """
    Purified E/B modes
    """

    npix = len(W)

    if Q.shape != (npix,) or U.shape != (npix,):
        raise ValueError("W, Q, U must all have shape (npix,)")

    nside = nside_from_npix(npix)

    ell = np.arange(lmax + 1, dtype=np.float64)

    n1 = np.zeros(lmax + 1, dtype=np.float64)
    n2 = np.zeros(lmax + 1, dtype=np.float64)

    n1[1:] = np.sqrt((ell[1:] + 1.0) * ell[1:])
    n2[2:] = np.sqrt((ell[2:] + 2.0) * (ell[2:]**2 - 1.0) * ell[2:])

    # Window alm
    wlm = sht.map2alm(lmax,lmax,W)

    # Compute del^1 W in spin-1 space
    tlm1_e = np.zeros_like(wlm)
    tlm1_b = np.zeros_like(wlm)
    for l in range(1, lmax + 1):
        tlm1_e[l,:] = wlm[l,:] * n1[l]
    W1 = sht.alm2map_spin(nside,1,np.array([tlm1_e, tlm1_b]), nthreads=0)

    # Compute del^2 W in spin-2 space
    tlm2_e = np.zeros_like(wlm)
    tlm2_b = np.zeros_like(wlm)
    for l in range(2, lmax + 1):
        tlm2_e[l,:] = wlm[l,:] * n2[l]
    W2 = sht.alm2map_spin(nside,2,np.array([tlm2_e, tlm2_b]), nthreads=0)

    # Inverse mask
    nonzero = W != 0.0

    W1_div = np.zeros_like(W1)
    W2_div = np.zeros_like(W2)
    W1_div[:,nonzero] = W1[:,nonzero] / W[None,nonzero]
    W2_div[:,nonzero] = W2[:,nonzero] / W[None,nonzero]

    W1 = W1_div
    W2 = W2_div

    # P2 = Q + iU
    P2_q = Q
    P2_u = U

    # P1 = conj(W1) * (Q + iU)
    P1_q = W1[0] * Q + W1[1] * U
    P1_u = W1[0] * U - W1[1] * Q

    # P0 = conj(W2) * (Q + iU)
    P0_q = W2[0] * Q + W2[1] * U
    P0_u = W2[0] * U - W2[1] * Q

    # Harmonic transforms
    alm2_e, alm2_b = sht.map2alm_spin(lmax,lmax,2,np.array([P2_q, P2_u]))
    alm1_e, alm1_b = sht.map2alm_spin(lmax,lmax,1,np.array([P1_q, P1_u]))

    # spin=0 case: transform each scalar component separately
    alm0_e = sht.map2alm(lmax,lmax,P0_q)
    alm0_b = sht.map2alm(lmax,lmax,P0_u)

    Elm = np.zeros_like(alm2_e, dtype=np.complex128)
    Blm = np.zeros_like(alm2_b, dtype=np.complex128)

    for l in range(2, lmax + 1):
        c1 = 2.0 * n1[l] / n2[l]
        c0 = 1.0 / n2[l]
        Elm[l,:] = alm2_e[l,:] + c1 * alm1_e[l,:] + c0 * alm0_e[l,:]
        Blm[l,:] = alm2_b[l,:] + c1 * alm1_b[l,:] + c0 * alm0_b[l,:]

    return Elm, Blm



def _angdist_from_chord(chord_distance):
    """
    Convert 3D chord distance on unit sphere to angular distance [rad].

    chord = 2 sin(theta / 2)
    theta = 2 arcsin(chord / 2)
    """
    x = np.clip(chord_distance / 2.0, 0.0, 1.0)
    return 2.0 * np.arcsin(x)


def _fill_small_holes_healpix(mask_nest, nside, hsize):
    """
    Fill small masked regions in NESTED ordering.

    Parameters
    ----------
    mask_nest : ndarray, shape (npix,)
        Integer/bool mask in NESTED order.
        1 = valid, 0 = masked.
    nside : int
        HEALPix nside.
    hsize : int
        Maximum number of pixels in a hole to fill.

    Returns
    -------
    filled : ndarray
        Mask after filling small holes.
    """
    if hsize <= 0:
        return mask_nest

    mask_nest = np.asarray(mask_nest, dtype=np.int8).copy()
    npix = mask_nest.size

    masked_pixels = np.flatnonzero(mask_nest == 0)
    visited = np.zeros(npix, dtype=bool)

    for start in masked_pixels:
        if visited[start] or mask_nest[start] != 0:
            continue

        # BFS over connected masked component
        component = []
        q = deque([start])
        visited[start] = True

        while q:
            p = q.popleft()
            component.append(p)

            neigh = hp.get_all_neighbours(nside, p, nest=True)
            for nb in neigh:
                if nb < 0:
                    continue
                if not visited[nb] and mask_nest[nb] == 0:
                    visited[nb] = True
                    q.append(nb)

            # Early exit: too large to fill
            if len(component) > hsize:
                # Still need to mark the whole connected component as visited.
                # Continue BFS, but we already know it will not be filled.
                pass

        if len(component) <= hsize:
            mask_nest[np.asarray(component, dtype=np.int64)] = 1

    return mask_nest


def apodize(rmask, ascale, order=1, holeminsize=0.0, fill_holes=True):
    """
    Parameters
    ----------
    rmask : ndarray, shape (npix,)
        Input window function.
        Pixels with rmask == 0 are treated as masked pixels.
    ascale : float
        Apodization length [deg] from the closest masked pixel.
    order : int, default 1
        Pixel ordering.
        1 means RING, otherwise NESTED, following the Fortran code.
    holeminsize : float, default 0.0
        Minimum hole size. Fortran comment says [arcmin], but the formula
        effectively treats this like an area compared with pixel area in arcmin^2.
    fill_holes : bool, default True
        If True, fill small holes before computing distance.

    Returns
    -------
    amask : ndarray, shape (npix,)
        Apodized window, returned in the same ordering as input.
    """
    rmask = np.asarray(rmask, dtype=np.float64)
    npix  = rmask.size
    nside = nside_from_npix(npix)

    if ascale <= 0:
        # No apodization length means keep original hard mask.
        return (rmask != 0.0).astype(np.float64) * rmask

    # Work internally in NESTED ordering
    if order == 1:
        rmask_nest = hp.reorder(rmask, r2n=True)
    else:
        rmask_nest = rmask.copy()

    # Integer mask: valid=1, masked=0
    mask_nest = (rmask_nest != 0.0).astype(np.int8)

    # Remove small holes
    if fill_holes and holeminsize > 0:
        # hsize = nint(holeminsize / pixel_area_arcmin2)
        pix_area_arcmin2 = (4.0 * np.pi / npix) * (60.0 * 180.0 / np.pi) ** 2
        hsize = int(np.rint(holeminsize / pix_area_arcmin2))

        if hsize > 0:
            mask_nest = _fill_small_holes_healpix(mask_nest, nside, hsize)

    # Compute angular distance from each valid pixel to nearest masked pixel.
    masked = np.flatnonzero(mask_nest == 0)
    valid = np.flatnonzero(mask_nest != 0)

    dist_nest = np.zeros(npix, dtype=np.float64)

    if masked.size == 0:
        # No masked pixels: distance is effectively infinite, window becomes 1.
        dist_nest[:] = np.inf
    elif valid.size > 0:
        # HEALPix pixel centers as 3D vectors.
        theta_m, phi_m = hp.pix2ang(nside, masked, nest=True)
        vec_m = np.column_stack(hp.ang2vec(theta_m, phi_m))

        tree = cKDTree(vec_m)

        theta_v, phi_v = hp.pix2ang(nside, valid, nest=True)
        vec_v = np.column_stack(hp.ang2vec(theta_v, phi_v))

        chord_dist, _ = tree.query(vec_v, k=1, workers=-1)
        dist_nest[valid] = _angdist_from_chord(chord_dist)

    # Back to original ordering before applying final formula
    if order == 1:
        dist = hp.reorder(dist_nest, n2r=True)
    else:
        dist = dist_nest

    # Compute apodization:
    #
    # x = (1 - cos(distance)) / (1 - cos(ascale))
    # y = min(1, sqrt(x))
    # amask = (y - sin(2*pi*y)/(2*pi)) * rmask
    denom = 1.0 - np.cos(np.deg2rad(ascale))

    x = (1.0 - np.cos(dist)) / denom
    y = np.minimum(1.0, np.sqrt(x))

    amask = (y - np.sin(2.0 * np.pi * y) / (2.0 * np.pi)) * rmask

    # Numerical cleanup
    amask = np.where(rmask == 0.0, 0.0, amask)

    return amask.astype(np.float64)


def map_mul_lfunc(imap, lfunc):
    """
    Convert map to alm, multiply alm by lfunc[l], then convert back to map.
    """
    npix = imap.size
    nside = nside_from_npix(npix)
    lmax = lfunc.size - 1

    alm = hp.map2alm(imap,lmax=lmax,mmax=lmax)
    alm = hp.almxfl(alm,lfunc[: lmax + 1],mmax=lmax,inplace=False)
    omap = hp.alm2map(alm,nside=nside,lmax=lmax,mmax=lmax,verbose=False)

    return np.asarray(omap, dtype=np.float64)


def cosin_healpix(nside: int, nest: bool = False):
    """
    Return cos(theta) as a function of HEALPix pixel index.

    Parameters
    ----------
    nside : int
        Number of HEALPix pixels, npix = 12*nside**2.
    nest : bool, default False
        If False, use RING ordering.
        If True, use NESTED ordering.

    Returns
    -------
    cosin : ndarray, shape (npix,)
        cos(theta) at each pixel center.
    """
    npix = 12*nside**2
    ipix = np.arange(npix, dtype=np.int64)
    theta, _ = hp.pix2ang(nside, ipix, nest=nest)

    return np.cos(theta).astype(np.float64)


def polcoord2angle(theta, phi, verbose=False):
    """
    The function is based on the fortran code provided by Takashi Hamana and Ryuichi Takahashi.

    Parameters
    ----------
    theta, phi : ndarray, shape (npix,)
        Source-plane coordinates [rad].
    verbose : bool
        Print diagnostic messages.

    Returns
    -------
    angle : ndarray, shape (2,npix)
        Deflection angle vector.
    """
    theta = np.asarray(theta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)

    if theta.shape != phi.shape:
        raise ValueError("theta and phi must have the same shape")

    npix  = theta.size
    nside = nside_from_npix(npix)

    if verbose:
        print("nside", nside)
        print("size:", theta.size, phi.size)
        print("obtain image plane theta/phi")

    ipix = np.arange(npix, dtype=np.int64)
    theta_i, phi_i = hp.pix2ang(nside, ipix, nest=False)

    if verbose:
        print("theta/phi to angle")

    deltaphi = phi - phi_i

    cosalp = np.cos(theta_i) * np.cos(theta) + np.sin(theta_i) * np.sin(theta) * np.cos(deltaphi)
    cosalp = np.clip(cosalp, -1.0, 1.0)

    alpha = np.arccos(cosalp)

    sinalp = np.sin(alpha)

    angle = np.zeros((2,npix), dtype=np.float64)

    nonzero = sinalp != 0.0

    if np.any(nonzero):
        ti = theta_i[nonzero]
        th = theta[nonzero]
        dp = deltaphi[nonzero]
        ca = cosalp[nonzero]
        sa = sinalp[nonzero]
        al = alpha[nonzero]

        denom = np.sin(ti) * sa

        good = denom != 0.0

        cosdelta = np.zeros_like(al)
        sindelta = np.zeros_like(al)

        cosdelta[good] = (np.cos(th[good]) - np.cos(ti[good]) * ca[good]) / denom[good]
        sindelta[good] = (np.sin(dp[good]) * np.sin(th[good])) / sa[good]

        idx = np.flatnonzero(nonzero)
        idx_good = idx[good]

        angle[0,idx_good] = -al[good] * cosdelta[good]
        angle[1,idx_good] =  al[good] * sindelta[good]

    return angle


def polcoord2angle_alm(lmax, theta, phi, verbose=False):
    """
    Parameters
    ----------
    lmax : int
        Maximum multipole.
    theta, phi : ndarray, shape (12*nside**2,)
        Source-plane coordinates [rad].
    verbose : bool
        Print diagnostic messages.

    Returns
    -------
    glm, clm : ndarray
        Gradient and curl modes.
    """
    theta = np.asarray(theta, dtype=np.float64)
    phi   = np.asarray(phi, dtype=np.float64)
    npix  = theta.size
    nside = nside_from_npix(npix)

    if verbose:
        print("convert to angle")

    angle = polcoord2angle(theta, phi, verbose=verbose)

    if verbose:
        print("deflection angle to its alms")

    # map2alm_spin expects [component_1, component_2]
    dlm_g, dlm_c = sht.map2alm_spin(lmax,lmax,1,angle)

    if verbose:
        print("convert to glm and clm")

    glm = np.zeros_like(dlm_g, dtype=np.complex128)
    clm = np.zeros_like(dlm_c, dtype=np.complex128)

    for l in range(1, lmax + 1):
        norm = np.sqrt(float(l * l + l))
        glm[l,:] = dlm_g[l,:] / norm
        clm[l,:] = dlm_c[l,:] / norm

    return glm, clm


def calc_mfs_from_derivatives(nu, der0, der1, der2):
    """
    Compute 2D Minkowski functionals from derivative maps.

    Parameters
    ----------
    nu : ndarray, shape (bn,)
        Threshold bins.
    der0 : ndarray, shape (npix,)
        Scalar field.
    der1 : ndarray, shape (npix, 2)
        First derivatives.
    der2 : ndarray, shape (npix, 3)
        Second derivatives.
        der2[:,0] = tt
        der2[:,1] = tp-like component
        der2[:,2] = pp-like component

    Returns
    -------
    V : ndarray, shape (bn, 3)
        V0, V1, V2.
    """
    npix = der0.size
    nside = nside_from_npix(npix)

    if nu.size < 2:
        raise ValueError("nu must contain at least two bins")

    dnu = nu[1] - nu[0]

    g1 = der1[0]
    g2 = der1[1]

    dtt = der2[0]
    dtp = der2[1]
    dpp = der2[2]

    grad2 = g1**2 + g2**2
    grad = np.sqrt(grad2)

    V = np.zeros((nu.size, 3), dtype=np.float64)

    for i, nui in enumerate(nu):
        # V0
        V[i, 0] = np.count_nonzero(der0 >= nui) / float(npix)
        shell = np.abs(der0 - nui) <= dnu / 2.0
        if not np.any(shell):
            continue

        # V1
        V[i, 1] = np.sum(0.25 / dnu * grad[shell]) / float(npix)

        # V2
        valid = shell & (grad2 != 0.0)

        if np.any(valid):
            numerator = 2.0 * g1[valid] * g2[valid] * dtp[valid] - g1[valid]**2 * dpp[valid] - g2[valid]**2 * dtt[valid]
            v2_pix = 0.5 / np.pi * (1.0 / dnu) * numerator / grad2[valid]
            V[i, 2] = np.sum(v2_pix) / float(npix)

    return V

    
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



