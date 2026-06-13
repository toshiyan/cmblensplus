from ._core import lib_rec_tau


def qtt(lmax, rlmin, rlmax, fC, Tlm1, Tlm2, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstruct inhomogeneous optical depth from the temperature quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output tau alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        TT spectrum, with bounds ``1:rlmax``.
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
    alm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Amplitude-modulation alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_tau.qtt(
        lmax, rlmin, rlmax, fC, Tlm1, Tlm2,
        nside_t, verbose, nthreads=nthreads
    )


def qeb(lmax, rlmin, rlmax, fCE, Elm, Blm, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstruct amplitude modulation from the EB quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output amplitude-modulation alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fCE : array_like of float, shape (rlmax + 1,)
        EE spectrum, with bounds ``0:rlmax``.
    Elm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered E-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Blm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered B-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    alm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Amplitude-modulation alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_tau.qeb(
        lmax, rlmin, rlmax, fCE, Elm, Blm,
        nside_t, verbose, nthreads=nthreads
    )
