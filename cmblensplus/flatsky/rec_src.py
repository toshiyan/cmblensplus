from . import libflatsky


def qtt(nx, ny, D, rL, T1, T2):
    """
    Reconstruct point-source fields from the temperature quadratic estimator.

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
    T1 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the first inverse-variance filtered
        temperature map.
    T2 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the second inverse-variance filtered
        temperature map.

    Returns
    -------
    slm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the point-source field.
    """
    return libflatsky.rec_src.qtt(nx, ny, D, rL, T1, T2)


def qte(nx, ny, D, rL, T, E):
    """
    Reconstruct point-source fields from the suboptimal TE quadratic estimator.

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
    T : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        temperature map.
    E : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        E-mode map.

    Returns
    -------
    slm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the point-source field.
    """
    return libflatsky.rec_src.qte(nx, ny, D, rL, T, E)


def qtb(nx, ny, D, rL, T, B):
    """
    Reconstruct point-source fields from the TB quadratic estimator.

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
    T : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        temperature map.
    B : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        B-mode map.

    Returns
    -------
    slm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the point-source field.
    """
    return libflatsky.rec_src.qtb(nx, ny, D, rL, T, B)


def qee(nx, ny, D, rL, E1, E2):
    """
    Reconstruct point-source fields from the EE quadratic estimator.

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
    E1 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the first inverse-variance filtered
        E-mode map.
    E2 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the second inverse-variance filtered
        E-mode map.

    Returns
    -------
    slm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the point-source field.
    """
    return libflatsky.rec_src.qee(nx, ny, D, rL, E1, E2)


def qeb(nx, ny, D, rL, E, B):
    """
    Reconstruct point-source fields from the EB quadratic estimator.

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
    E : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        E-mode map.
    B : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        B-mode map.

    Returns
    -------
    slm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the point-source field.
    """
    return libflatsky.rec_src.qeb(nx, ny, D, rL, E, B)