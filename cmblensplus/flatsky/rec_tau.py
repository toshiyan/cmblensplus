from . import libflatsky


def qtt(nx, ny, D, rL, fC, T1, T2):
    """
    Reconstruct patchy tau from the temperature quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        Temperature power spectrum on the 2D Fourier grid.
    T1 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the first inverse-variance filtered
        temperature field.
    T2 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the second inverse-variance filtered
        temperature field.

    Returns
    -------
    tlm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of patchy tau.
    """
    return libflatsky.rec_tau.qtt(nx, ny, D, rL, fC, T1, T2)


def qte(nx, ny, D, rL, fC, T, E):
    """
    Reconstruct patchy tau from the suboptimal TE quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        TE cross-power spectrum on the 2D Fourier grid.
    T : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        temperature field.
    E : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        E-mode field.

    Returns
    -------
    tlm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of patchy tau.
    """
    return libflatsky.rec_tau.qte(nx, ny, D, rL, fC, T, E)


def qtb(nx, ny, D, rL, fC, T, B):
    """
    Reconstruct patchy tau from the TB quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        TE cross-power spectrum on the 2D Fourier grid.
    T : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        temperature field.
    B : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        B-mode field.

    Returns
    -------
    tlm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of patchy tau.
    """
    return libflatsky.rec_tau.qtb(nx, ny, D, rL, fC, T, B)


def qee(nx, ny, D, rL, fC, E1, E2):
    """
    Reconstruct patchy tau from the EE quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fC : ndarray of float, shape (nx, ny)
        EE power spectrum on the 2D Fourier grid.
    E1 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the first inverse-variance filtered
        E-mode field.
    E2 : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the second inverse-variance filtered
        E-mode field.

    Returns
    -------
    tlm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of patchy tau.
    """
    return libflatsky.rec_tau.qee(nx, ny, D, rL, fC, E1, E2)


def qeb(nx, ny, D, rL, fE, fB, E, B):
    """
    Reconstruct patchy tau from the EB quadratic estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    fE : ndarray of float, shape (nx, ny)
        EE power spectrum on the 2D Fourier grid.
    fB : ndarray of float, shape (nx, ny)
        BB power spectrum on the 2D Fourier grid.
    E : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        E-mode field.
    B : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of the inverse-variance filtered
        B-mode field.

    Returns
    -------
    tlm : ndarray of complex, shape (nx, ny)
        Two-dimensional Fourier modes of patchy tau.
    """
    return libflatsky.rec_tau.qeb(nx, ny, D, rL, fE, fB, E, B)