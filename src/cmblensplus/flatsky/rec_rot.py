from . import libflatsky


def qte(nx, ny, D, rL, fC, T, E):
    """
    Reconstruct anisotropic polarization-rotation angles from the TE quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalently ``dLx / (2 * pi)`` and
        ``dLy / (2 * pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        TE cross-power spectrum on the 2D grid.
    T : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the inverse-variance filtered temperature map.
    E : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the inverse-variance filtered E-mode map.

    Returns
    -------
    alm : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the anisotropic polarization-rotation angles.
    """
    return libflatsky.rec_rot.qte(nx, ny, D, rL, fC, T, E)


def qtb(nx, ny, D, rL, fC, T, B):
    """
    Reconstruct anisotropic polarization-rotation angles from the TB quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalently ``dLx / (2 * pi)`` and
        ``dLy / (2 * pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        TE cross-power spectrum on the 2D grid.
    T : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the inverse-variance filtered temperature map.
    B : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the inverse-variance filtered B-mode map.

    Returns
    -------
    alm : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the anisotropic polarization-rotation angles.
    """
    return libflatsky.rec_rot.qtb(nx, ny, D, rL, fC, T, B)


def qee(nx, ny, D, rL, fC, E1, E2):
    """
    Reconstruct anisotropic polarization-rotation angles from the EE quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalently ``dLx / (2 * pi)`` and
        ``dLy / (2 * pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        EE power spectrum on the 2D grid.
    E1 : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the first inverse-variance filtered E-mode map.
    E2 : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the second inverse-variance filtered E-mode map.

    Returns
    -------
    alm : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the anisotropic polarization-rotation angles.
    """
    return libflatsky.rec_rot.qee(nx, ny, D, rL, fC, E1, E2)


def qeb(nx, ny, D, rL, EE, E, B, BB=0):
    """
    Reconstruct anisotropic polarization-rotation angles from the EB quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalently ``dLx / (2 * pi)`` and
        ``dLy / (2 * pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    EE : ndarray of float, shape (nx, ny)
        EE power spectrum on the 2D grid.
    E : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the inverse-variance filtered E-mode map.
    B : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the inverse-variance filtered B-mode map.
    BB : ndarray of float or float, optional
        Theory B-mode spectrum on the 2D grid. If a scalar is given, it is
        applied as a constant. Default is 0.

    Returns
    -------
    alm : ndarray of complex, shape (nx, ny)
        2D Fourier modes of the anisotropic polarization-rotation angles.
    """
    return libflatsky.rec_rot.qeb(nx, ny, D, rL, EE, E, B, BB)
    