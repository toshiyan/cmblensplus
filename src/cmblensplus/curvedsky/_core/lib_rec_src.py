from . import lib_utils as utils
from . import sht_ducc0 as sht
import numpy as np


def qtt(lmax, rlmin, rlmax, Tlm1, Tlm2, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstructing point sources from the temperature quadratic estimator.
    """

    Tlm1 = utils._validate_alm("Tlm1", Tlm1, rlmax)
    Tlm2 = utils._validate_alm("Tlm2", Tlm2, rlmax)

    nside = utils._default_nside(lmax, nside_t)
    if verbose:
        print(f"calc src-TT estimator with nside= {nside}")

    # Tlm1 -> scalar map
    alm1 = utils._fill_component("scal", 0, Tlm1, rlmin, rlmax)
    map1 = sht.alm2map(nside, alm1[0], nthreads=nthreads)

    # Tlm2 -> scalar map
    alm2 = utils._fill_component("scal", 0, Tlm2, rlmin, rlmax)
    map2 = sht.alm2map(nside, alm2[0], nthreads=nthreads)

    # convolution:
    map0 = map1 * map2

    # scalar map -> alm
    slm = sht.map2alm(lmax, lmax, map0, nthreads=nthreads)

    return 0.5*slm
