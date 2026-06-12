from . import lib_utils as utils
from . import sht_ducc0 as sht
import numpy as np


def qtt(lmax, rlmin, rlmax, fC, Tlm1, Tlm2, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstructing amplitude modulation from the temperature quadratic estimator.
    """

    fC   = utils._validate_cl("fC", fC, rlmax)
    Tlm1 = utils._validate_alm("Tlm1", Tlm1, rlmax)
    Tlm2 = utils._validate_alm("Tlm2", Tlm2, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc tau-TT estimator with nside= {nside}")

    # Tlm1 -> scalar map
    alm1 = utils._fill_component("scal", 0, Tlm1, rlmin, rlmax)
    map1 = sht.alm2map(nside, alm1[0], nthreads=nthreads)

    # fC * Tlm2 -> scalar map
    alm2 = utils._fill_component("scal", 0, Tlm2, rlmin, rlmax, fC)
    map2 = sht.alm2map(nside, alm2[0], nthreads=nthreads)

    # convolution:
    map0 = map1 * map2
    tlm = sht.map2alm(lmax, lmax, map0, nthreads=nthreads)

    return tlm


def qeb(lmax, rlmin, rlmax, fCE, Elm, Blm, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstructing amplitude modulation from the EB quadratic estimator.
    """

    fCE = utils._validate_cl("fCE", fCE, rlmax)
    Elm = utils._validate_alm("Elm", Elm, rlmax)
    Blm = utils._validate_alm("Blm", Blm, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc tau-EB estimator with nside= {nside}")

    # B spin-2 map
    almB = utils._fill_component("spin", 1, Blm, rlmin, rlmax)
    A = sht.alm2map_spin(nside, 2, almB, nthreads=nthreads)

    # fCE * E spin-2 map
    almE = utils._fill_component("spin", 0, Elm, rlmin, rlmax, fCE)
    A2 = sht.alm2map_spin(nside, 2, almE, nthreads=nthreads)

    # convolution:
    map0 = A[0] * A2[0] + A[1] * A2[1]
    tlm = sht.map2alm(lmax, lmax, map0, nthreads=nthreads)

    return tlm
