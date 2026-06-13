from . import libbasic


def lintemplate(lmax, elmin, elmax, klmin, klmax, CE, Cm, WE, Wm, gtype='p'):
    """
    Estimate the lensing-template B-mode power spectrum.

    Wiener filters are given as inputs.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output spectrum.
    elmin : int
        Minimum multipole of the E-mode spectrum.
    elmax : int
        Maximum multipole of the E-mode spectrum.
    klmin : int
        Minimum multipole of the lensing-potential or convergence spectrum.
    klmax : int
        Maximum multipole of the lensing-potential or convergence spectrum.
    CE : array_like of float, shape (dlmax + 1,)
        E-mode power spectrum, with bounds ``0:dlmax``.
    Cm : array_like of float, shape (dlmax + 1,)
        Lensing-potential or convergence power spectrum, with bounds
        ``0:dlmax``.
    WE : array_like of float, shape (dlmax + 1,)
        Wiener filter for the E-mode, with bounds ``0:dlmax``.
    Wm : array_like of float, shape (dlmax + 1,)
        Wiener filter for the lensing potential or convergence, with bounds
        ``0:dlmax``.
    gtype : {'p', 'k'}, optional
        Type of input ``Cm`` and ``Wm``. Use ``'p'`` for lensing potential
        or ``'k'`` for convergence. Default is ``'p'``.

    Returns
    -------
    CB : ndarray of float, shape (lmax + 1,)
        Lensing B-mode power spectrum, with bounds ``0:lmax``.
    """
    return libbasic.delens.lintemplate(
        lmax, elmin, elmax, klmin, klmax, CE, Cm, WE, Wm, gtype
    )


def lensingbb(lmax, dlmin, dlmax, CE, Cp):
    """
    Compute the lensing B-mode power spectrum as a convolution of the E-mode
    and lensing-potential power spectra.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output spectrum.
    dlmin : int
        Minimum multipole of the E-mode and lensing-potential spectra used for
        delensing.
    dlmax : int
        Maximum multipole of the E-mode and lensing-potential spectra used for
        delensing.
    CE : array_like of float, shape (dlmax + 1,)
        E-mode power spectrum, with bounds ``0:dlmax``.
    Cp : array_like of float, shape (dlmax + 1,)
        Lensing-potential power spectrum, with bounds ``0:dlmax``.

    Returns
    -------
    CB : ndarray of float, shape (lmax + 1,)
        Lensing B-mode power spectrum, with bounds ``0:lmax``.
    """
    return libbasic.delens.lensingbb(lmax, dlmin, dlmax, CE, Cp)


def delensbias_dom(lmax, dlmin, dlmax, CE, CB, Cp, NP1, NP2, Ag):
    """
    Compute the dominant term of the delensing bias in B-mode internal
    delensing.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output spectrum.
    dlmin : int
        Minimum multipole of the E-mode and lensing-potential spectra used for
        delensing.
    dlmax : int
        Maximum multipole of the E-mode and lensing-potential spectra used for
        delensing.
    CE : array_like of float, shape (dlmax + 1,)
        E-mode power spectrum, with bounds ``0:dlmax``.
    CB : array_like of float, shape (dlmax + 1,)
        B-mode power spectrum, with bounds ``0:dlmax``.
    Cp : array_like of float, shape (dlmax + 1,)
        Lensing-potential power spectrum, with bounds ``0:dlmax``.
    NP1 : array_like of float, shape (dlmax + 1,)
        Polarization noise spectrum for lensing reconstruction, with bounds
        ``0:dlmax``.
    NP2 : array_like of float, shape (dlmax + 1,)
        Polarization noise spectrum for the B-mode to be delensed, with bounds
        ``0:dlmax``.
    Ag : array_like of float, shape (dlmax + 1,)
        Lensing-reconstruction noise, with bounds ``0:dlmax``.

    Returns
    -------
    DB : ndarray of float, shape (lmax + 1,)
        Dominant delensing-bias contribution to the B-mode power spectrum,
        with bounds ``0:lmax``.
    """
    return libbasic.delens.delensbias_dom(
        lmax, dlmin, dlmax, CE, CB, Cp, NP1, NP2, Ag
    )
    