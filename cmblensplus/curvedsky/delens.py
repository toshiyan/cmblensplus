from ._core import lib_delens


def lensingb(lmax, elmin, elmax, plmin, plmax, wElm, wplm, nside_t=0, gtype='p'):
    """
    Compute the lensing B-mode alm as a convolution of the Wiener-filtered
    E-mode alm and the lensing-potential alm.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing B-mode alm.
    elmin : int
        Minimum multipole of the Wiener-filtered E-mode alm.
    elmax : int
        Maximum multipole of the Wiener-filtered E-mode alm.
    plmin : int
        Minimum multipole of the Wiener-filtered lensing-potential alm.
    plmax : int
        Maximum multipole of the Wiener-filtered lensing-potential alm.
    wElm : ndarray of complex, shape (elmax + 1, elmax + 1)
        Wiener-filtered E-mode alm, with bounds ``(0:elmax, 0:elmax)``.
    wplm : ndarray of complex, shape (plmax + 1, plmax + 1)
        Wiener-filtered lensing-potential or convergence alm, with bounds
        ``(0:plmax, 0:plmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : {'p', 'k'}, optional
        Type of input ``wplm``. Use ``'p'`` for lensing potential
        :math:`\phi`, or ``'k'`` for convergence :math:`\kappa`.
        Default is ``'p'``.

    Returns
    -------
    lBlm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Lensing B-mode alm, with bounds ``(0:lmax, 0:lmax)``.
    """
    return lib_delens.lensingb(
        lmax, elmin, elmax, plmin, plmax, wElm, wplm, nside_t, gtype
    )


def shiftvec(nside, plm, nremap):
    r"""
    Return the anti-deflection vector :math:`\beta` at each Healpix pixel.

    The anti-deflection vector is used for delensing and satisfies

    .. math::

        \beta(\hat{n}) + \alpha^{\rm iw}(\hat{n} + \beta(\hat{n})) = 0,

    where :math:`\alpha^{\rm iw}` is the filtered lensing deflection vector.
    See arXiv:1701.01712.

    Parameters
    ----------
    nside : int
        Nside of the output shift-vector map.
    plm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Wiener-filtered lensing-potential alm, with bounds
        ``(0:lmax, 0:lmax)``.
    nremap : int
        Number of iterations used to compute the shift vector.

    Returns
    -------
    beta : ndarray of float, shape (12 * nside**2, 2)
        Two-dimensional shift vector, with bounds ``(0:npix-1, 1:2)``.
    """
    return lib_delens.shiftvec(nside, plm, nremap)


def phi2grad(nside, plm):
    """
    Return the deflection vector at each Healpix pixel.

    Parameters
    ----------
    nside : int
        Nside of the output deflection-vector map.
    plm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.

    Returns
    -------
    grad : ndarray of float, shape (12 * nside**2, 2)
        Two-dimensional deflection vector, with bounds ``(0:npix-1, 1:2)``.
    """
    return lib_delens.phi2grad(nside, plm)
