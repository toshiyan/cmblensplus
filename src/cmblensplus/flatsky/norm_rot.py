from . import libflatsky


def qeb(nx, ny, D, rL, IE, IB, EE, eL, BB=0):
    """
    Return the normalization of the EB quadratic estimator for anisotropic
    polarization-rotation angles.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, or equivalently ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    IE : ndarray of float, shape (nx, ny)
        Inverse of the observed E-mode power spectrum on the 2D grid.
    IB : ndarray of float, shape (nx, ny)
        Inverse of the observed B-mode power spectrum on the 2D grid.
    EE : ndarray of float, shape (nx, ny)
        Theoretical E-mode power spectrum on the 2D grid.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.
    BB : float or ndarray of float, optional
        Theoretical B-mode power spectrum on the 2D grid. If an array is
        given, its shape should be ``(nx, ny)``. Default is 0.

    Returns
    -------
    Aa : ndarray of complex, shape (nx, ny)
        Normalization of anisotropic polarization-rotation angles on the 2D
        grid.
    """
    return libflatsky.norm_rot.qeb(nx, ny, D, rL, IE, IB, EE, eL, BB)
    