from ._core import lib_rec_src


def qtt(lmax, rlmin, rlmax, Tlm1, Tlm2, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstruct point sources from the temperature quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output point-source alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    Tlm1 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        First inverse-variance filtered temperature alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Tlm2 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Second inverse-variance filtered temperature alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    slm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Point-source alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_src.qtt(
        lmax, rlmin, rlmax, Tlm1, Tlm2,
        nside_t, verbose, nthreads=nthreads
    )
