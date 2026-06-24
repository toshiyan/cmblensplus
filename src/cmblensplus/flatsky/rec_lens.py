from . import libflatsky


def qtt(nx, ny, D, rL, fC, T1, T2, gtype=''):
    """
    Reconstruct the CMB lensing potential and curl mode from the temperature
    quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        Temperature power spectrum on the 2D Fourier grid.
    T1 : ndarray of complex, shape (nx, ny)
        Fourier modes of the first inverse-variance filtered temperature map.
    T2 : ndarray of complex, shape (nx, ny)
        Fourier modes of the second inverse-variance filtered temperature map.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    glm : ndarray of complex, shape (nx, ny)
        Fourier modes of the CMB lensing potential.
    clm : ndarray of complex, shape (nx, ny)
        Fourier modes of the curl mode, or pseudo lensing potential.
    """
    return libflatsky.rec_lens.qtt(nx, ny, D, rL, fC, T1, T2, gtype)


def qte(nx, ny, D, rL, fC, T, E, gtype=''):
    """
    Reconstruct the CMB lensing potential and curl mode from the suboptimal TE
    quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        TE cross-power spectrum on the 2D Fourier grid.
    T : ndarray of complex, shape (nx, ny)
        Fourier modes of the inverse-variance filtered temperature map.
    E : ndarray of complex, shape (nx, ny)
        Fourier modes of the inverse-variance filtered E-mode map.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    glm : ndarray of complex, shape (nx, ny)
        Fourier modes of the CMB lensing potential.
    clm : ndarray of complex, shape (nx, ny)
        Fourier modes of the curl mode, or pseudo lensing potential.
    """
    return libflatsky.rec_lens.qte(nx, ny, D, rL, fC, T, E, gtype)


def qtb(nx, ny, D, rL, fC, T, B, gtype=''):
    """
    Reconstruct the CMB lensing potential and curl mode from the TB quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        TE cross-power spectrum on the 2D Fourier grid.
    T : ndarray of complex, shape (nx, ny)
        Fourier modes of the inverse-variance filtered temperature map.
    B : ndarray of complex, shape (nx, ny)
        Fourier modes of the inverse-variance filtered B-mode map.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    glm : ndarray of complex, shape (nx, ny)
        Fourier modes of the CMB lensing potential.
    clm : ndarray of complex, shape (nx, ny)
        Fourier modes of the curl mode, or pseudo lensing potential.
    """
    return libflatsky.rec_lens.qtb(nx, ny, D, rL, fC, T, B, gtype)


def qee(nx, ny, D, rL, fC, E1, E2, gtype=''):
    """
    Reconstruct the CMB lensing potential and curl mode from the EE quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        EE power spectrum on the 2D Fourier grid.
    E1 : ndarray of complex, shape (nx, ny)
        Fourier modes of the first inverse-variance filtered E-mode map.
    E2 : ndarray of complex, shape (nx, ny)
        Fourier modes of the second inverse-variance filtered E-mode map.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    glm : ndarray of complex, shape (nx, ny)
        Fourier modes of the CMB lensing potential.
    clm : ndarray of complex, shape (nx, ny)
        Fourier modes of the curl mode, or pseudo lensing potential.
    """
    return libflatsky.rec_lens.qee(nx, ny, D, rL, fC, E1, E2, gtype)


def qeb(nx, ny, D, rL, fC, E, B, gtype=''):
    """
    Reconstruct the CMB lensing potential and curl mode from the EB quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        EE power spectrum on the 2D Fourier grid.
    E : ndarray of complex, shape (nx, ny)
        Fourier modes of the inverse-variance filtered E-mode map.
    B : ndarray of complex, shape (nx, ny)
        Fourier modes of the inverse-variance filtered B-mode map.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    glm : ndarray of complex, shape (nx, ny)
        Fourier modes of the CMB lensing potential.
    clm : ndarray of complex, shape (nx, ny)
        Fourier modes of the curl mode, or pseudo lensing potential.
    """
    return libflatsky.rec_lens.qeb(nx, ny, D, rL, fC, E, B, gtype)


def qbb(nx, ny, D, rL, fC, B1, B2, gtype=''):
    """
    Reconstruct the CMB lensing potential and curl mode from the BB quadratic
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points in the x direction.
    ny : int
        Number of Fourier grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        BB power spectrum on the 2D Fourier grid.
    B1 : ndarray of complex, shape (nx, ny)
        Fourier modes of the first inverse-variance filtered B-mode map.
    B2 : ndarray of complex, shape (nx, ny)
        Fourier modes of the second inverse-variance filtered B-mode map.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    glm : ndarray of complex, shape (nx, ny)
        Fourier modes of the CMB lensing potential.
    clm : ndarray of complex, shape (nx, ny)
        Fourier modes of the curl mode, or pseudo lensing potential.
    """
    return libflatsky.rec_lens.qbb(nx, ny, D, rL, fC, B1, B2, gtype)