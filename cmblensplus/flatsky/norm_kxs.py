from . import libflatsky


def qtt(nx, ny, D, rL, OT, TT, eL):
    """
    Return the normalization of the temperature quadratic estimators for the
    CMB lensing potential and patchy optical depth.

    Parameters
    ----------
    nx : int
        Number of grid points along the x direction.
    ny : int
        Number of grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2 * pi)`` and
        ``dLy / (2 * pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    OT : ndarray of float, shape (nx, ny)
        Inverse of the observed temperature power spectrum on the 2D grid,
        with bounds ``(0:nx-1, 0:ny-1)``.
    TT : ndarray of float, shape (nx, ny)
        Theoretical temperature power spectrum on the 2D grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Aks : ndarray of complex, shape (nx, ny)
        Lensing-tau cross-normalization on the 2D grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    """
    return libflatsky.norm_kxs.qtt(nx, ny, D, rL, OT, TT, eL)
    