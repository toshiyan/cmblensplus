from . import libflatsky


def map2alm(nx, ny, D, map):
    """
    Perform a discrete Fourier transform of a 2D map.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    map : ndarray of float, shape (nx, ny)
        Input map on a 2D grid.

    Returns
    -------
    alm : ndarray of complex, shape (nx, ny)
        Fourier modes on the 2D grid.
    """
    return libflatsky.utils.map2alm(nx, ny, D, map)


def alm2map(nx, ny, D, alm):
    """
    Perform an inverse discrete Fourier transform of 2D Fourier modes.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    alm : ndarray of complex, shape (nx, ny)
        Fourier modes on a 2D grid.

    Returns
    -------
    map : ndarray of float, shape (nx, ny)
        Transformed map on the 2D grid.
    """
    return libflatsky.utils.alm2map(nx, ny, D, alm)


def el2d(nx, ny, D):
    r"""
    Return the absolute value of the multipole on a 2D Fourier grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.

    Returns
    -------
    els : ndarray of float, shape (nx, ny)
        Absolute value of the Fourier mode,
        ``sqrt(Lx**2 + Ly**2)``.
    """
    return libflatsky.utils.el2d(nx, ny, D)


def elarrays(nx, ny, D):
    r"""
    Return the Fourier-grid coordinates and their radial quantities.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.

    Returns
    -------
    elx : ndarray of float, shape (nx, ny)
        Fourier coordinate ``Lx``.
    ely : ndarray of float, shape (nx, ny)
        Fourier coordinate ``Ly``.
    els : ndarray of float, shape (nx, ny)
        Absolute value of the Fourier mode,
        ``sqrt(Lx**2 + Ly**2)``.
    eli : ndarray of float, shape (nx, ny)
        Inverse of ``els``.
    """
    return libflatsky.utils.elarrays(nx, ny, D)


def elmask(nx, ny, D, lmin=0, lmax=1000, lxcut=0, lycut=0):
    """
    Return a mask in 2D Fourier space.

    The mask is unity where ``lmin <= |L| <= lmax``,
    ``|Lx| >= lxcut``, and ``|Ly| >= lycut``. It is zero elsewhere.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    lmin : int, optional
        Minimum multipole. Default is 0.
    lmax : int, optional
        Maximum multipole. Default is 1000.
    lxcut : int, optional
        Remove modes with ``|Lx| < lxcut``. Default is 0.
    lycut : int, optional
        Remove modes with ``|Ly| < lycut``. Default is 0.

    Returns
    -------
    lmask : ndarray of float, shape (nx, ny)
        Fourier-space mask.
    """
    return libflatsky.utils.elmask(nx, ny, D, lmin, lmax, lxcut, lycut)


def ulm_flat(nx, ny, D, ulm, lmax=0, alpha=1.0, deapod=True):
    """
    Convert spherical harmonic coefficients to a 2D flat-sky grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    ulm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Harmonic coefficients.
    lmax : int, optional
        Maximum multipole of ``ulm``. If 0, it is inferred from ``ulm``.
        Default is 0.
    alpha : float, optional
        Scaling parameter passed to the backend routine. Default is 1.0.
    deapod : bool, optional
        Whether to apply deapodization. Default is True.

    Returns
    -------
    ul2d : ndarray of complex, shape (nx, ny)
        Harmonic coefficients converted onto the 2D grid.
    """
    if lmax == 0:
        lmax = len(ulm[:, 0]) - 1
    return libflatsky.utils.ulm_flat(nx, ny, D, lmax, ulm, alpha, deapod)


def ulm2ulphi(ulm, lmax=0):
    """
    Convert spherical harmonic coefficients to ``u(l, phi_l)``.

    Parameters
    ----------
    ulm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Harmonic coefficients.
    lmax : int, optional
        Maximum multipole of ``ulm``. If 0, it is inferred from ``ulm``.
        Default is 0.

    Returns
    -------
    ulphi : ndarray of complex, shape (lmax + 1, lmax + 1)
        Converted harmonic coefficients on the ``(l, phi_l)`` grid.
    """
    if lmax == 0:
        lmax = len(ulm[:, 0]) - 1
    return libflatsky.utils.ulm2ulphi(lmax, ulm)


def lphi_to_cartesian(nx, ny, D, lmax, ulphi):
    """
    Convert coefficients from an ``(l, phi_l)`` grid to a Cartesian 2D grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths.
    lmax : int
        Maximum multipole of ``ulphi``.
    ulphi : ndarray of complex
        Input coefficients on the ``(l, phi_l)`` grid.

    Returns
    -------
    array_like
        Coefficients converted to a Cartesian 2D grid.
    """
    return libflatsky.utils.lphi_to_cartesian(nx, ny, D, lmax, ulphi)


def deapodization_bilinear(nx, ny, D, imap, alpha):
    """
    Apply bilinear deapodization to a 2D map.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths.
    imap : ndarray of float, shape (nx, ny)
        Input map.
    alpha : float
        Deapodization parameter.

    Returns
    -------
    omap : ndarray of float, shape (nx, ny)
        Deapodized output map.
    """
    return libflatsky.utils.deapodization_bilinear(nx, ny, D, imap, alpha)


def alm2bcl(bn, oL, nx, ny, D, alm1, alm2=None, spc=''):
    """
    Compute a binned angular power spectrum from Fourier modes.

    Parameters
    ----------
    bn : int
        Number of multipole bins.
    oL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output spectrum.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    alm1 : ndarray of complex, shape (nx, ny)
        First Fourier-mode array.
    alm2 : ndarray of complex, shape (nx, ny), optional
        Second Fourier-mode array. If not given, ``alm1`` is used.
    spc : str, optional
        Multipole bin spacing. Use ``''`` for linear spacing or ``'log'``
        for logarithmic spacing. Default is ``''``.

    Returns
    -------
    Cb : ndarray of float, shape (bn,)
        Binned angular power spectrum.
    """
    if alm2 is None:
        alm2 = alm1
    return libflatsky.utils.alm2bcl(bn, oL, nx, ny, D, alm1, alm2, spc)


def c2d2bcl(nx, ny, D, c2d, bn, oL, spc=''):
    """
    Compute a binned 1D angular power spectrum from a 2D power spectrum.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    c2d : ndarray of float, shape (nx, ny)
        Input 2D power spectrum.
    bn : int
        Number of multipole bins.
    oL : array_like of int, shape (2,)
        Minimum and maximum multipoles of the output spectrum.
    spc : str, optional
        Multipole bin spacing. Use ``''`` for linear spacing or ``'log'``
        for logarithmic spacing. Default is ``''``.

    Returns
    -------
    Cb : ndarray of float, shape (bn,)
        Binned angular power spectrum.
    """
    return libflatsky.utils.c2d2bcl(nx, ny, D, c2d, bn, oL, spc)


def cl2c2d(nx, ny, D, lmin, lmax, Cl, method='linear'):
    """
    Interpolate a 1D angular power spectrum onto a 2D Fourier grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    lmin : int
        Minimum multipole of ``Cl`` to be interpolated.
    lmax : int
        Maximum multipole of ``Cl`` to be interpolated.
    Cl : array_like of float, shape (lmax + 1,)
        Input 1D angular power spectrum.
    method : str, optional
        Interpolation method. Use ``'linear'`` for linear interpolation or
        ``'step'`` for step interpolation. Default is ``'linear'``.

    Returns
    -------
    c2d : ndarray of float, shape (nx, ny)
        Interpolated 2D power spectrum.
    """
    return libflatsky.utils.cl2c2d(nx, ny, D, lmin, lmax, Cl, method)


def cb2c2d(bn, bc, nx, ny, D, lmin, lmax, Cb, method=''):
    """
    Interpolate a binned 1D angular power spectrum onto a 2D Fourier grid.

    Parameters
    ----------
    bn : int
        Number of multipole bins.
    bc : array_like of float, shape (bn,)
        Multipole bin centers.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    lmin : int
        Minimum multipole to be interpolated.
    lmax : int
        Maximum multipole to be interpolated.
    Cb : array_like of float, shape (bn,)
        Binned 1D angular power spectrum.
    method : str, optional
        Interpolation method from binned to unbinned spectrum. Use ``''``
        for spline interpolation or ``'linear'`` for linear interpolation.
        Default is ``''``.

    Returns
    -------
    c2d : ndarray of float, shape (nx, ny)
        Interpolated 2D power spectrum.
    """
    return libflatsky.utils.cb2c2d(bn, bc, nx, ny, D, lmin, lmax, Cb, method)


def gauss1alm(nx, ny, D, lmin, lmax, Cl):
    """
    Generate a random Gaussian field in 2D Fourier space.

    The generated field satisfies the Hermitian condition
    ``a_l.conjugate() = a_{-l}``.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    lmin : int
        Minimum multipole of ``Cl`` to be interpolated.
    lmax : int
        Maximum multipole of ``Cl`` to be interpolated.
    Cl : array_like of float, shape (lmax + 1,)
        Input 1D angular power spectrum.

    Returns
    -------
    alm : ndarray of complex, shape (nx, ny)
        Random Gaussian field on the 2D Fourier plane.
    """
    return libflatsky.utils.gauss1alm(nx, ny, D, lmin, lmax, Cl)


def gauss2alm(nx, ny, D, lmin, lmax, TT, TE, EE):
    """
    Generate two correlated random Gaussian fields in 2D Fourier space.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    lmin : int
        Minimum multipole of the spectra to be interpolated.
    lmax : int
        Maximum multipole of the spectra to be interpolated.
    TT : array_like of float, shape (lmax + 1,)
        First 1D auto power spectrum.
    TE : array_like of float, shape (lmax + 1,)
        1D cross power spectrum.
    EE : array_like of float, shape (lmax + 1,)
        Second 1D auto power spectrum.

    Returns
    -------
    tlm : ndarray of complex, shape (nx, ny)
        First random Gaussian field on the 2D Fourier plane.
    elm : ndarray of complex, shape (nx, ny)
        Second random Gaussian field on the 2D Fourier plane.
    """
    return libflatsky.utils.gauss2alm(nx, ny, D, lmin, lmax, TT, TE, EE)


def window_sin(nx, ny, D, ap=1, cut=1):
    """
    Return a sine window function.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    ap : float, optional
        Apodization parameter. The apodized range is
        ``(1 - ap) * cut * mapsize``. It ranges from 0, full apodization,
        to 1, no apodization. Default is 1.
    cut : float, optional
        Map cut parameter. The cut map size is ``cut * mapsize``.
        It ranges from 0, full cut, to 1, no cut. Default is 1.

    Returns
    -------
    W : ndarray of float, shape (nx, ny)
        Window function.
    """
    return libflatsky.utils.window_sin(nx, ny, D, ap, cut)


def window_norm(nx, ny, wind, num):
    """
    Return the normalization factor of a window function.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    wind : ndarray of float, shape (nx, ny)
        Window function.
    num : int
        Power of the window used for the normalization.

    Returns
    -------
    norm : float
        Window normalization factor.
    """
    return libflatsky.utils.window_norm(nx, ny, wind, num)


def window_norm_x(nx, ny, W1, W2, num):
    """
    Return the cross-normalization factor of two window functions.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    W1 : ndarray of float, shape (nx, ny)
        First window function.
    W2 : ndarray of float, shape (nx, ny)
        Second window function.
    num : int
        Power used for the normalization.

    Returns
    -------
    norm : float
        Cross-window normalization factor.
    """
    return libflatsky.utils.window_norm_x(nx, ny, W1, W2, num)


def rotation(nx, ny, rot, QU, rtype='f'):
    """
    Rotate a two-component field on a 2D grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    rot : ndarray of float, shape (nx, ny)
        Rotation angle map.
    QU : ndarray of float
        Input two-component field.
    rtype : str
        Rotation type ``l'' (small angle limit) or ``f'' (no assumption)

    Returns
    -------
    rotated : ndarray
        Rotated two-component field.
    """
    return libflatsky.utils.rotation(nx, ny, rot, QU, rtype)


def get_angle(nx, ny, D):
    """
    Return the Fourier-space angle on a 2D grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,) 
        Map side lengths.

    Returns
    -------
    angle : ndarray of float, shape (nx, ny)
        Fourier-space angle on the 2D grid.
    """
    return libflatsky.utils.get_angle(nx, ny, D)


def cutmap(ox, oy, cx, cy, omap):
    """
    Cut out a rectangular region from a 2D map.

    Parameters
    ----------
    ox : int
        Number of grid points in the x direction of the original map.
    oy : int
        Number of grid points in the y direction of the original map.
    cx : int
        Number of grid points in the x direction of the cut map.
    cy : int
        Number of grid points in the y direction of the cut map.
    omap : ndarray
        Original input map.

    Returns
    -------
    cmap : ndarray
        Cut-out map.
    """
    return libflatsky.utils.cutmap(ox, oy, cx, cy, omap)
    