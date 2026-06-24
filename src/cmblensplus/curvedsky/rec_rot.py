from ._core import lib_rec_rot


def qeb(lmax, rlmin, rlmax, fCE, Elm, Blm, nside_t=0, verbose=False, nthreads=0):
    """
    Reconstruct the polarization rotation angle from the EB quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output rotation-angle alm.
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
        Rotation-angle alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_rot.qeb(
        lmax, rlmin, rlmax, fCE, Elm, Blm,
        nside_t, verbose, nthreads=nthreads
    )
