from . import libflatsky


def qtt(nx, ny, D, rL, OT, TT, eL):
    """
    Return the normalization of the temperature quadratic estimator for the
    CMB lensing potential and curl mode.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    OT : ndarray of float, shape (nx, ny)
        Inverse observed temperature power spectrum on the 2D grid, with
        bounds ``(nx, ny)``.
    TT : ndarray of float, shape (nx, ny)
        Theoretical temperature power spectrum on the 2D grid, with bounds
        ``(nx, ny)``.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.qtt(nx, ny, D, rL, OT, TT, eL)


def n0tt(nx, ny, D, rL, OT0, OT1, TT, eL):
    """
    Return the disconnected noise normalization of the temperature quadratic
    estimator for the CMB lensing potential and curl mode.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    OT0 : ndarray of float, shape (nx, ny)
        First inverse observed temperature power spectrum on the 2D grid.
    OT1 : ndarray of float, shape (nx, ny)
        Second inverse observed temperature power spectrum on the 2D grid.
    TT : ndarray of float, shape (nx, ny)
        Theoretical temperature power spectrum on the 2D grid.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.n0tt(nx, ny, D, rL, OT0, OT1, TT, eL)


def n0ttc(nx, ny, D, rL, OT0, OT1, TT, eL):
    """
    Return the disconnected noise normalization of the temperature quadratic
    estimator for the CMB lensing potential and curl mode, allowing complex
    inverse observed spectra.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    OT0 : ndarray of complex, shape (nx, ny)
        First inverse observed temperature power spectrum on the 2D grid.
    OT1 : ndarray of complex, shape (nx, ny)
        Second inverse observed temperature power spectrum on the 2D grid.
    TT : ndarray of float, shape (nx, ny)
        Theoretical temperature power spectrum on the 2D grid.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.n0ttc(nx, ny, D, rL, OT0, OT1, TT, eL)


def qte(nx, ny, D, rL, OT, OE, TE, eL):
    """
    Return the normalization of the TE quadratic estimator for the CMB lensing
    potential and curl mode.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    OT : ndarray of float, shape (nx, ny)
        Inverse observed temperature power spectrum on the 2D grid.
    OE : ndarray of float, shape (nx, ny)
        Inverse observed E-mode power spectrum on the 2D grid.
    TE : ndarray of float, shape (nx, ny)
        Theoretical TE cross spectrum on the 2D grid.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.qte(nx, ny, D, rL, OT, OE, TE, eL)


def qtb(nx, ny, D, OT, OB, TE, rL, eL):
    """
    Return the normalization of the TB quadratic estimator for the CMB lensing
    potential and curl mode.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    OT : ndarray of float, shape (nx, ny)
        Inverse observed temperature power spectrum on the 2D grid.
    OB : ndarray of float, shape (nx, ny)
        Inverse observed B-mode power spectrum on the 2D grid.
    TE : ndarray of float, shape (nx, ny)
        Theoretical TE cross spectrum on the 2D grid.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.qtb(nx, ny, D, OT, OB, TE, rL, eL)


def qee(nx, ny, D, OE, EE, rL, eL):
    """
    Return the normalization of the EE quadratic estimator for the CMB lensing
    potential and curl mode.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    OE : ndarray of float, shape (nx, ny)
        Inverse observed E-mode power spectrum on the 2D grid.
    EE : ndarray of float, shape (nx, ny)
        Theoretical E-mode power spectrum on the 2D grid.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.qee(nx, ny, D, OE, EE, rL, eL)


def qeb(nx, ny, D, OE, OB, EE, rL, eL):
    """
    Return the normalization of the EB quadratic estimator for the CMB lensing
    potential and curl mode.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``, with bounds ``0:1``.
    OE : ndarray of float, shape (nx, ny)
        Inverse observed E-mode power spectrum on the 2D grid.
    OB : ndarray of float, shape (nx, ny)
        Inverse observed B-mode power spectrum on the 2D grid.
    EE : ndarray of float, shape (nx, ny)
        Theoretical E-mode power spectrum on the 2D grid.
    rL : array_like of int, shape (2,)
        Minimum and maximum CMB multipoles used for reconstruction.
    eL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output normalization spectrum.

    Returns
    -------
    Ag : ndarray of complex, shape (nx, ny)
        Normalization of the CMB lensing potential on the 2D grid.
    Ac : ndarray of complex, shape (nx, ny)
        Normalization of the curl mode, or pseudo lensing potential, on the
        2D grid.
    """
    return libflatsky.norm_lens.qeb(nx, ny, D, OE, OB, EE, rL, eL)
    