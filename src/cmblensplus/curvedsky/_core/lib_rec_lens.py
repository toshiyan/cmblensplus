from . import lib_utils as utils
from . import sht_ducc0 as sht
import numpy as np


def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t=0,gtype='',verbose=False,nthreads=0):
    """
    Reconstructing CMB lensing potential and its curl mode from the temperature quadratic estimator
    """

    fC   = utils._validate_cl("fC", fC, rlmax)
    Tlm1 = utils._validate_alm("Tlm1", Tlm1, rlmax)
    Tlm2 = utils._validate_alm("Tlm2", Tlm2, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qTT lens estimator with nside= {nside}")

    ilk = utils._ilk(lmax, gtype)

    alm1 = utils._fill_component('scal', 0, Tlm1, rlmin, rlmax)
    amap = sht.alm2map(nside, alm1[0], nthreads=nthreads)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    weights = fC * np.sqrt(ell * (ell + 1.0))
    alm1 = utils._fill_component('spin', 0, Tlm2, rlmin, rlmax, weights)
    bmap = sht.alm2map_spin(nside, 1, alm1, nthreads=nthreads)
    
    # convolution
    bmap[0] *= amap
    bmap[1] *= amap
    blm = sht.map2alm_spin(lmax, lmax, 1, bmap, nthreads=nthreads)

    # compute glm and clm
    return utils._finish_gradient_curl(lmax, blm, ilk, prefactor=1.0)


def qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside_t=0,gtype='',verbose=False,nthreads=0):
    """
    Reconstructing CMB lensing potential and its curl mode from the TE quadratic estimator
    """
    fC  = utils._validate_cl("fC", fC, rlmax)
    Tlm = utils._validate_alm("Tlm", Tlm, rlmax)
    Elm = utils._validate_alm("Elm", Elm, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qTE lens estimator with nside= {nside}")

    ilk = utils._ilk(lmax, gtype)

    # first part
    alm1 = utils._fill_component('spin', 0, Elm, rlmin, rlmax)
    A = sht.alm2map_spin(nside, 2, alm1, nthreads=nthreads)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    x1 = (ell + 2.0) * (ell - 1.0)
    x3 = (ell - 2.0) * (ell + 3.0)

    s1 = np.zeros_like(x1, dtype=float)
    s3 = np.zeros_like(x3, dtype=float)

    np.sqrt(x1, out=s1, where=(x1 >= 0.0))
    np.sqrt(x3, out=s3, where=(x3 >= 0.0))

    w1 = fC * s1
    w3 = fC * s3

    alm1 = utils._fill_component('spin', 0, Tlm, rlmin, rlmax, w1)
    alm3 = utils._fill_component('spin', 0, Tlm, rlmin, rlmax, w3)
    A1 = sht.alm2map_spin(nside, 1, alm1, nthreads=nthreads)
    A3 = sht.alm2map_spin(nside, 3, alm3, nthreads=nthreads)

    # second part
    alm0 = utils._fill_component('scal', 0, Tlm, rlmin, rlmax)
    AT = sht.alm2map(nside, alm0[0], nthreads=nthreads)

    wAE = fC * np.sqrt(ell * (ell + 1.0))
    elm1 = utils._fill_component('spin', 0, Elm, rlmin, rlmax, wAE)
    AE = sht.alm2map_spin(nside, 1, elm1, nthreads=nthreads)

    # convolution
    maps = utils._common_spin_product(A, A1, A3)
    maps[0] += 2.0 * AT * AE[0]
    maps[1] += 2.0 * AT * AE[1]
    blm = sht.map2alm_spin(lmax, lmax, 1, maps, nthreads=nthreads)

    return utils._finish_gradient_curl(lmax, blm, ilk, prefactor=0.5)

    
def qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside_t=0,gtype='',verbose=False,nthreads=0):
    """
    Reconstructing CMB lensing potential and its curl mode from the TB quadratic estimator
    """
    fC  = utils._validate_cl("fC", fC, rlmax)
    Tlm = utils._validate_alm("Tlm", Tlm, rlmax)
    Blm = utils._validate_alm("Blm", Blm, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qTB lens estimator with nside= {nside}")

    ilk = utils._ilk(lmax, gtype)

    alm1 = utils._fill_component('spin', 1, Blm, rlmin, rlmax)
    A = sht.alm2map_spin(nside, 2, alm1, nthreads=nthreads)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    x1 = (ell + 2.0) * (ell - 1.0)
    x3 = (ell - 2.0) * (ell + 3.0)

    s1 = np.zeros_like(x1, dtype=float)
    s3 = np.zeros_like(x3, dtype=float)

    np.sqrt(x1, out=s1, where=(x1 >= 0.0))
    np.sqrt(x3, out=s3, where=(x3 >= 0.0))

    w1 = fC * s1
    w3 = fC * s3
    alm1 = utils._fill_component('spin', 0, Tlm, rlmin, rlmax, w1)
    alm3 = utils._fill_component('spin', 0, Tlm, rlmin, rlmax, w3)
    A1 = sht.alm2map_spin(nside, 1, alm1, nthreads=nthreads)
    A3 = sht.alm2map_spin(nside, 3, alm3, nthreads=nthreads)

    # convolution
    maps = utils._common_spin_product(A, A1, A3)
    zlm  = sht.map2alm_spin(lmax, lmax, 1, maps, nthreads=nthreads)

    return utils._finish_gradient_curl(lmax, zlm, ilk, prefactor=0.5)

    
def qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside_t=0,gtype='',verbose=False,nthreads=0):
    """
    Reconstructing CMB lensing potential and its curl mode from the EE quadratic estimator
    """
    fC   = utils._validate_cl("fC", fC, rlmax)
    Elm1 = utils._validate_alm("Elm1", Elm1, rlmax)
    Elm2 = utils._validate_alm("Elm2", Elm2, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qEE lens estimator with nside= {nside}")

    ilk = utils._ilk(lmax, gtype)

    alm = utils._fill_component('spin', 0, Elm1, rlmin, rlmax)
    A = sht.alm2map_spin(nside, 2, alm, nthreads=nthreads)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    x1 = (ell + 2.0) * (ell - 1.0)
    x3 = (ell - 2.0) * (ell + 3.0)

    s1 = np.zeros_like(x1, dtype=float)
    s3 = np.zeros_like(x3, dtype=float)

    np.sqrt(x1, out=s1, where=(x1 >= 0.0))
    np.sqrt(x3, out=s3, where=(x3 >= 0.0))

    w1 = fC * s1
    w3 = fC * s3

    alm = utils._fill_component('spin', 0, Elm2, rlmin, rlmax, w1)
    blm = utils._fill_component('spin', 0, Elm2, rlmin, rlmax, w3)
    A1  = sht.alm2map_spin(nside, 1, alm, nthreads=nthreads)
    A3  = sht.alm2map_spin(nside, 3, blm, nthreads=nthreads)

    # convolution
    maps = utils._common_spin_product(A, A1, A3)
    alm_out = sht.map2alm_spin(lmax, lmax, 1, maps, nthreads=nthreads)

    return utils._finish_gradient_curl(lmax, alm_out, ilk, prefactor=0.5)


def qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside_t=0,gtype='',verbose=False,nthreads=0):
    """
    Reconstructing CMB lensing potential and its curl mode from the EB quadratic estimator
    """
    fC  = utils._validate_cl("fC", fC, rlmax)
    Elm = utils._validate_alm("Elm", Elm, rlmax)
    Blm = utils._validate_alm("Blm", Blm, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qEB lens estimator with nside= {nside}")

    ilk = utils._ilk(lmax, gtype)

    alm1 = utils._fill_component('spin', 1, Blm, rlmin, rlmax)
    A = sht.alm2map_spin(nside, 2, alm1, nthreads=nthreads)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    x1 = (ell + 2.0) * (ell - 1.0)
    x3 = (ell - 2.0) * (ell + 3.0)

    s1 = np.zeros_like(x1, dtype=float)
    s3 = np.zeros_like(x3, dtype=float)

    np.sqrt(x1, out=s1, where=(x1 >= 0.0))
    np.sqrt(x3, out=s3, where=(x3 >= 0.0))

    w1 = fC * s1
    w3 = fC * s3
    alm1 = utils._fill_component('spin', 0, Elm, rlmin, rlmax, w1)
    alm3 = utils._fill_component('spin', 0, Elm, rlmin, rlmax, w3)
    A1 = sht.alm2map_spin(nside, 1, alm1, nthreads=nthreads)
    A3 = sht.alm2map_spin(nside, 3, alm3, nthreads=nthreads)

    # convolution
    maps = utils._common_spin_product(A, A1, A3)
    tlm = sht.map2alm_spin(lmax, lmax, 1, maps, nthreads=nthreads)

    return utils._finish_gradient_curl(lmax, tlm, ilk, prefactor=0.5)

    
def qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside_t=0,gtype='',verbose=False,nthreads=0):
    """
    Reconstructing CMB lensing potential and its curl mode from the BB quadratic estimator
    """
    fC   = utils._validate_cl("fC", fC, rlmax)
    Blm1 = utils._validate_alm("Blm1", Blm1, rlmax)
    Blm2 = utils._validate_alm("Blm2", Blm2, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc qBB lens estimator with nside= {nside}")

    ilk = utils._ilk(lmax, gtype)

    alm1 = utils._fill_component('spin', 1, Blm1, rlmin, rlmax)
    A = sht.alm2map_spin(nside, 2, alm1, nthreads=nthreads)

    ell = np.arange(rlmax + 1, dtype=np.float64)
    x1 = (ell + 2.0) * (ell - 1.0)
    x3 = (ell - 2.0) * (ell + 3.0)

    s1 = np.zeros_like(x1, dtype=float)
    s3 = np.zeros_like(x3, dtype=float)

    np.sqrt(x1, out=s1, where=(x1 >= 0.0))
    np.sqrt(x3, out=s3, where=(x3 >= 0.0))

    w1 = fC * s1
    w3 = fC * s3
    alm1 = utils._fill_component('spin', 1, Blm2, rlmin, rlmax, w1)
    alm3 = utils._fill_component('spin', 1, Blm2, rlmin, rlmax, w3)
    A1 = sht.alm2map_spin(nside, 1, alm1, nthreads=nthreads)
    A3 = sht.alm2map_spin(nside, 3, alm3, nthreads=nthreads)

    # convolution
    maps = utils._common_spin_product(A, A1, A3)
    tlm = sht.map2alm_spin(lmax, lmax, 1, maps, nthreads=nthreads)

    return utils._finish_gradient_curl(lmax, tlm, ilk, prefactor=0.5)


__all__ = [
    "qtt",
    "qte",
    "qtb",
    "qee",
    "qeb",
    "qbb",
]
