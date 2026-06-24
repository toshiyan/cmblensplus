from . import libflatsky


def qtt(nx, ny, D, rL, OT, TT, eL):
    """
    Return the normalization of the temperature quadratic estimator between
    the CMB lensing potential and patchy optical depth.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2 * pi)`` and
        ``dLy / (2 * pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    OT : ndarray of float, shape (nx, ny)
        Inverse observed temperature power spectrum on the 2D grid.
    TT : ndarray of float, shape (nx, ny)
        Theoretical temperature power spectrum on the 2D grid.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Akt : ndarray of complex, shape (nx, ny)
        Lensing-tau cross-normalization on the 2D grid.
    """
    return libflatsky.norm_kxt.qtt(nx, ny, D, rL, OT, TT, eL)
    