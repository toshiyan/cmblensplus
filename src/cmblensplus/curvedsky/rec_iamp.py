from ._core import lib_rec_iamp


def qeb(lmax, rlmin, rlmax, EB, Elm, Blm, nside_t=0, verbose=False):
    """
    Reconstruct amplitude modulation with the odd EB quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output amplitude-modulation alm.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    EB : array_like of float, shape (rlmax + 1,)
        EB spectrum, with bounds ``0:rlmax``.
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

    Returns
    -------
    alm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Reconstructed amplitude-modulation alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_iamp.qeb(lmax, rlmin, rlmax, EB, Elm, Blm, nside_t, verbose)
