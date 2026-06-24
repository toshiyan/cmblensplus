from . import lib_utils as utils
from . import sht_ducc0 as sht
import numpy as np


def qeb(lmax, rlmin, rlmax, fC, Elm, Blm, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstructing CMB polarization rotation angle from the EB quadratic estimator.
    """

    fC  = utils._validate_cl("fC", fC, rlmax)
    Elm = utils._validate_alm("Elm", Elm, rlmax)
    Blm = utils._validate_alm("Blm", Blm, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc rot-EB estimator with nside= {nside}")

    # B spin-2 map
    almB = utils._fill_component("spin", 1, Blm, rlmin, rlmax)
    AB = sht.alm2map_spin(nside, 2, almB, nthreads=nthreads)

    # fC * E spin-2 map
    almE = utils._fill_component("spin", 0, Elm, rlmin, rlmax, fC)
    AE = sht.alm2map_spin(nside, 2, almE, nthreads=nthreads)

    # convolution:
    map0 = AB[0] * AE[1] - AB[1] * AE[0]

    # scalar map -> alm
    alm = sht.map2alm(lmax, lmax, map0, nthreads=nthreads)

    return -2*alm

__all__ = [
    #"qte",
    #"qtb",
    #"qee",
    "qeb",
    #"qbb",
]
