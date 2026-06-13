from . import libcurvedsky as lib_norm
# to be replaced with pytempura


def qte(est, lmax, rlmin, rlmax, TB, OCT, OCE, lfac=''):
    """
    Compute the normalization of the reconstructed imaginary CMB lensing
    potential and curl mode from the TE quadratic estimator.

    Parameters
    ----------
    est : str
        Estimator type.
    lmax : int
        Maximum multipole of the output normalization spectrum.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    TB : array_like of float, shape (rlmax + 1,)
        Theory TB spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Output type. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_imag.qte(est, lmax, rlmin, rlmax, TB, OCT, OCE, lfac)


def qtb(est, lmax, rlmin, rlmax, TB, OCT, OCB, lfac=''):
    """
    Compute the normalization of the reconstructed imaginary CMB lensing
    potential and curl mode from the TB quadratic estimator.

    Parameters
    ----------
    est : str
        Estimator type.
    lmax : int
        Maximum multipole of the output normalization spectrum.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    TB : array_like of float, shape (rlmax + 1,)
        Theory TE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Output type. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_imag.qtb(est, lmax, rlmin, rlmax, TB, OCT, OCB, lfac)


def qee(est, lmax, rlmin, rlmax, fC, OCE, lfac=''):
    """
    Compute the normalization of the reconstructed imaginary CMB lensing
    potential and curl mode from the E-mode quadratic estimator.

    Parameters
    ----------
    est : str
        Estimator type.
    lmax : int
        Maximum multipole of the output normalization spectrum.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Output type. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_imag.qee(est, lmax, rlmin, rlmax, fC, OCE, lfac)


def qeb(est, lmax, rlmin, rlmax, fC, OCE, OCB, lfac=''):
    """
    Compute the normalization of the reconstructed imaginary CMB lensing
    potential and curl mode from the EB quadratic estimator.

    Parameters
    ----------
    est : str
        Estimator type.
    lmax : int
        Maximum multipole of the output normalization spectrum.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        Theory EB spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Output type. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_imag.qeb(est, lmax, rlmin, rlmax, fC, OCE, OCB, lfac)


def qbb(est, lmax, rlmin, rlmax, fC, OCB, lfac=''):
    """
    Compute the normalization of the reconstructed imaginary CMB lensing
    potential and curl mode from the B-mode quadratic estimator.

    Parameters
    ----------
    est : str
        Estimator type.
    lmax : int
        Maximum multipole of the output normalization spectrum.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        Theory BB spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Output type. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_imag.qbb(est, lmax, rlmin, rlmax, fC, OCB, lfac)


def qbb_asym(est, lmax, rlmin, rlmax, fC, OCB1, OCB2, lfac=''):
    """
    Compute the asymmetric normalization of the reconstructed imaginary CMB
    lensing potential and curl mode from the B-mode quadratic estimator.

    Parameters
    ----------
    est : str
        Estimator type.
    lmax : int
        Maximum multipole of the output normalization spectrum.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        Theory BB spectrum, with bounds ``0:rlmax``.
    OCB1 : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum for experiment 1, with bounds ``0:rlmax``.
    OCB2 : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum for experiment 2, with bounds ``0:rlmax``.
    lfac : str, optional
        Output type. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_imag.qbb_asym(
        est, lmax, rlmin, rlmax, fC, OCB1, OCB2, lfac
    )