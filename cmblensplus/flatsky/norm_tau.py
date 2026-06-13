from . import libflatsky
import numpy


def qtt(nx, ny, D, rL, OT, TT, eL):
    """
    Return the normalization of the temperature quadratic estimator for
    patchy optical depth.

    Parameters
    ----------
    nx : int
        Number of grids along the x direction.
    ny : int
        Number of grids along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
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
    At : ndarray of complex, shape (nx, ny)
        Normalization of patchy optical depth on the 2D grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    """
    return libflatsky.norm_tau.qtt(nx, ny, D, rL, OT, TT, eL)


def qeb(nx, ny, D, rL, IE, IB, EE, eL):
    """
    Return the normalization of the EB quadratic estimator for patchy optical
    depth.

    Parameters
    ----------
    nx : int
        Number of grids along the x direction.
    ny : int
        Number of grids along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    IE : ndarray of float, shape (nx, ny)
        Inverse of the observed E-mode power spectrum on the 2D grid, with
        bounds ``(0:nx-1, 0:ny-1)``.
    IB : ndarray of float, shape (nx, ny)
        Inverse of the observed B-mode power spectrum on the 2D grid, with
        bounds ``(0:nx-1, 0:ny-1)``.
    EE : ndarray of float, shape (nx, ny)
        Theoretical E-mode power spectrum on the 2D grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    At : ndarray of complex, shape (nx, ny)
        Normalization of patchy optical depth on the 2D grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    """
    return libflatsky.norm_tau.qeb(nx, ny, D, rL, IE, IB, EE, eL)