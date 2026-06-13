# This module provides the normalization of the quadratic estimators

import numpy as np
from . import libcurvedsky as lib_norm
# to be replaced with pytempura


def qtt(est, lmax, rlmin, rlmax, TT, OCT, lfac=''):
    """
    Compute the normalization of reconstructed fields from the temperature
    quadratic estimator.

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
    TT : array_like of float, shape (rlmax + 1,)
        Theory TT spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalizations. The first component is for the lensing-potential
        estimator and the second component is for the curl-mode estimator.
    """
    return lib_norm.norm_quad.qtt(est, lmax, rlmin, rlmax, TT, OCT, lfac)


def qte(est, lmax, rlmin, rlmax, TE, OCT, OCE, lfac=''):
    """
    Compute the normalization of reconstructed fields from the TE quadratic
    estimator.

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
    TE : array_like of float, shape (rlmax + 1,)
        Theory TE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalizations. The first component is for the lensing-potential
        estimator and the second component is for the curl-mode estimator.
    """
    return lib_norm.norm_quad.qte(est, lmax, rlmin, rlmax, TE, OCT, OCE, lfac)


def qtb(est, lmax, rlmin, rlmax, TE, OCT, OCB, lfac=''):
    """
    Compute the normalization of reconstructed fields from the TB quadratic
    estimator.

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
    TE : array_like of float, shape (rlmax + 1,)
        Theory TE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalizations. The first component is for the lensing-potential
        estimator and the second component is for the curl-mode estimator.
    """
    return lib_norm.norm_quad.qtb(est, lmax, rlmin, rlmax, TE, OCT, OCB, lfac)


def qee(est, lmax, rlmin, rlmax, EE, OCE, lfac=''):
    """
    Compute the normalization of reconstructed fields from the EE quadratic
    estimator.

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
    EE : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_quad.qee(est, lmax, rlmin, rlmax, EE, OCE, lfac)


def qeb(est, lmax, rlmin, rlmax, EE, OCE, OCB, BB=0, lfac=''):
    """
    Compute the normalization of reconstructed fields from the EB quadratic
    estimator.

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
    EE : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    BB : array_like of float, shape (rlmax + 1,), optional
        Theory BB spectrum, with bounds ``0:rlmax``. If set to 0, a zero
        spectrum with the same shape as ``EE`` is used.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    if BB == 0:
        BB = 0. * EE
    return lib_norm.norm_quad.qeb(est, lmax, rlmin, rlmax, EE, OCE, OCB, BB, lfac)


def qbb(est, lmax, rlmin, rlmax, BB, OCB, lfac=''):
    """
    Compute the normalization of reconstructed fields from the BB quadratic
    estimator.

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
    BB : array_like of float, shape (rlmax + 1,)
        Theory BB spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Al : ndarray of float, shape (2, lmax + 1)
        Normalization, with bounds ``(0:1, 0:lmax)``.
    """
    return lib_norm.norm_quad.qbb(est, lmax, rlmin, rlmax, BB, OCB, lfac)


def qttte(est, lmax, rlmin, rlmax, fCTT, fCTE, OCT, OCE, OCTE, lfac=''):
    """
    Compute the correlation between unnormalized TT and TE quadratic
    estimators.

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
    fCTT : array_like of float, shape (rlmax + 1,)
        Theory TT spectrum, with bounds ``0:rlmax``.
    fCTE : array_like of float, shape (rlmax + 1,)
        Theory TE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCTE : array_like of float, shape (rlmax + 1,)
        Observed TE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Ig : ndarray of float, shape (lmax + 1,)
        Correlation between lensing-potential estimators, with bounds
        ``0:lmax``.
    Ic : ndarray of float, shape (lmax + 1,)
        Correlation between curl-mode estimators, with bounds ``0:lmax``.
    """
    return lib_norm.norm_quad.qttte(
        est, lmax, rlmin, rlmax, fCTT, fCTE, OCT, OCE, OCTE, lfac
    )


def qttee(est, lmax, rlmin, rlmax, fCTT, fCEE, OCT, OCE, OCTE, lfac=''):
    """
    Compute the correlation between unnormalized TT and EE quadratic
    estimators.

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
    fCTT : array_like of float, shape (rlmax + 1,)
        Theory TT spectrum, with bounds ``0:rlmax``.
    fCEE : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCTE : array_like of float, shape (rlmax + 1,)
        Observed TE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Ig : ndarray of float, shape (lmax + 1,)
        Correlation between lensing-potential estimators, with bounds
        ``0:lmax``.
    Ic : ndarray of float, shape (lmax + 1,)
        Correlation between curl-mode estimators, with bounds ``0:lmax``.
    """
    return lib_norm.norm_quad.qttee(
        est, lmax, rlmin, rlmax, fCTT, fCEE, OCT, OCE, OCTE, lfac
    )


def qteee(est, lmax, rlmin, rlmax, fCEE, fCTE, OCT, OCE, OCTE, lfac=''):
    """
    Compute the correlation between unnormalized TE and EE quadratic
    estimators.

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
    fCEE : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    fCTE : array_like of float, shape (rlmax + 1,)
        Theory TE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCTE : array_like of float, shape (rlmax + 1,)
        Observed TE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Ig : ndarray of float, shape (lmax + 1,)
        Correlation between lensing-potential estimators, with bounds
        ``0:lmax``.
    Ic : ndarray of float, shape (lmax + 1,)
        Correlation between curl-mode estimators, with bounds ``0:lmax``.
    """
    return lib_norm.norm_quad.qteee(
        est, lmax, rlmin, rlmax, fCEE, fCTE, OCT, OCE, OCTE, lfac
    )


def qtbeb(est, lmax, rlmin, rlmax, fCEE, fCBB, fCTE, OCT, OCE, OCB, OCTE, lfac=''):
    """
    Compute the correlation between unnormalized TB and EB quadratic
    estimators.

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
    fCEE : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    fCBB : array_like of float, shape (rlmax + 1,)
        Theory BB spectrum, with bounds ``0:rlmax``.
    fCTE : array_like of float, shape (rlmax + 1,)
        Theory TE spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    OCTE : array_like of float, shape (rlmax + 1,)
        Observed TE spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Ig : ndarray of float, shape (lmax + 1,)
        Correlation between lensing-potential estimators, with bounds
        ``0:lmax``.
    Ic : ndarray of float, shape (lmax + 1,)
        Correlation between curl-mode estimators, with bounds ``0:lmax``.
    """
    return lib_norm.norm_quad.qtbeb(
        est, lmax, rlmin, rlmax, fCEE, fCBB, fCTE, OCT, OCE, OCB, OCTE, lfac
    )


def qmv(lmax, QDO, Al, Il):
    """
    Compute the minimum-variance estimator normalization.

    Currently, the BB estimator is ignored.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output power spectra.
    QDO : array_like of bool, shape (6,)
        Flags specifying which estimators are combined for the
        minimum-variance estimator. The order is TT, TE, EE, TB, EB, and BB.
        Currently, BB is always False.
    Al : array_like of float, shape (5, lmax + 1)
        Normalizations of each estimator: TT, TE, EE, TB, and EB.
    Il : array_like of float, shape (4, lmax + 1)
        Correlations between different estimators: TT x TE, TT x EE,
        TE x EE, and TB x EB.

    Returns
    -------
    MV : ndarray of float, shape (lmax + 1,)
        Normalization of the minimum-variance estimator, with bounds
        ``0:lmax``.
    Nl : ndarray of float, shape (6, lmax + 1)
        Weights for each estimator: TT, TE, EE, TB, EB, and BB, with bounds
        ``(0:5, 0:lmax)``. The BB weight is zero.
    """
    return lib_norm.norm_quad.qmv(lmax, QDO, Al, Il)


def qall(est, QDO, lmax, rlmin, rlmax, fC, OC, lfac=''):
    """
    Compute all quadratic-estimator normalizations and MV weights.

    Currently, the BB estimator is ignored.

    Parameters
    ----------
    est : str
        Estimator type.
    QDO : array_like of bool, shape (6,)
        Flags specifying which estimators are combined for the
        minimum-variance estimator. The order is TT, TE, EE, TB, EB, and BB.
    lmax : int
        Maximum multipole of the output power spectra.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (4, rlmax + 1)
        Theory CMB angular power spectra, with components TT, EE, BB, and TE,
        and bounds ``(0:3, 0:rlmax)``.
    OC : array_like of float, shape (4, rlmax + 1)
        Observed CMB angular power spectra, with components TT, EE, BB, and
        TE, and bounds ``(0:3, 0:rlmax)``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Ag : ndarray of float, shape (6, lmax + 1)
        Normalizations of the TT, TE, EE, TB, EB, and MV estimators for the
        lensing potential, with bounds ``(0:5, 0:lmax)``.
    Ac : ndarray of float, shape (6, lmax + 1)
        Same as ``Ag``, but for the curl mode.
    Nlg : ndarray of float, shape (6, lmax + 1)
        Weights for TT, TE, EE, TB, EB, and BB estimators for the lensing
        potential, with bounds ``(0:5, 0:lmax)``. The BB weight is zero.
    Nlc : ndarray of float, shape (6, lmax + 1)
        Same as ``Nlg``, but for the curl mode.
    """
    return lib_norm.norm_quad.qall(est, QDO, lmax, rlmin, rlmax, fC, OC, lfac)


def qeb_iter(lmax, elmax, rlmin, rlmax, dlmin, dlmax, CE, OCE, OCB, Cpp, iter=1, conv=0.00001):
    """
    Compute the iterative EB normalization of the reconstructed CMB lensing
    potential and curl mode.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output normalization.
    elmax : int
        Maximum multipole of the input EE spectra ``CE`` and ``OCE``.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    dlmin : int
        Minimum multipole of the E mode and lensing potential used for
        delensing.
    dlmax : int
        Maximum multipole of the E mode and lensing potential used for
        delensing.
    CE : array_like of float, shape (elmax + 1,)
        Theory EE angular power spectrum, with bounds ``0:elmax``.
    OCE : array_like of float, shape (elmax + 1,)
        Observed EE spectrum, with bounds ``0:elmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    Cpp : array_like of float, shape (dlmax + 1,)
        Theory lensing-potential spectrum, with bounds ``0:dlmax``.
    iter : int, optional
        Number of iterations. Default is 1, corresponding to no iteration.
    conv : float, optional
        Convergence threshold for the iteration. Default is 0.00001.

    Returns
    -------
    Ag : ndarray of float, shape (lmax + 1,)
        CMB lensing-potential normalization, with bounds ``0:lmax``.
    Ac : ndarray of float, shape (lmax + 1,)
        Curl-mode, or pseudo lensing-potential, normalization, with bounds
        ``0:lmax``.
    """
    return lib_norm.norm_quad.qeb_iter(
        lmax, elmax, rlmin, rlmax, dlmin, dlmax, CE, OCE, OCB, Cpp, iter, conv
    )


def qall_iter(lmax, rlmin, rlmax, ucl, lcl, ocl, eb='qe'):
    """
    Compute reconstruction-noise curves for several estimator combinations.

    The returned rows correspond to TT, TE, EE, TB, EB, EE+EB, and
    TT+TE+EE+EB.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output reconstruction-noise curves.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    ucl : array_like of float
        Unlensed CMB angular power spectra.
    lcl : array_like of float
        Lensed CMB angular power spectra.
    ocl : array_like of float
        Observed CMB angular power spectra.
    eb : {'qe', 'iter'}, optional
        EB reconstruction method. Use ``'qe'`` for the quadratic estimator,
        or any other value for the iterative EB reconstruction. Default is
        ``'qe'``.

    Returns
    -------
    Ag : ndarray of float, shape (7, lmax + 1)
        Lensing-potential reconstruction-noise curves.
    Ac : ndarray of float, shape (7, lmax + 1)
        Curl-mode reconstruction-noise curves.
    """
    Ag = np.zeros((7, lmax + 1))
    Ac = np.zeros((7, lmax + 1))

    # QDO = TT+TE+EE, and Ag[5] = TT+TE+EE
    Ag[:6, :], Ac[:6, :], nlg, nlc = qall(
        'lens', [True, True, True, False, False, False],
        lmax, rlmin, rlmax, lcl, ocl
    )
    Ag[3, :], Ac[3, :] = qtb(
        'lens', lmax, rlmin, rlmax, lcl[3, :], ocl[0, :], ocl[2, :]
    )
    if eb == 'qe':
        # QE reconstruction for EB
        Ag[4, :], Ac[4, :] = qeb(
            'lens', lmax, rlmin, rlmax, lcl[1, :], ocl[1, :], ocl[2, :]
        )
    else:
        # Iterative reconstruction for EB
        Ag[4, :], Ac[4, :] = qeb_iter(
            lmax, rlmax, rlmin, rlmax, rlmin, rlmax,
            ucl[1, :], ocl[1, :], ocl[2, :], ucl[3, :],
            iter=10, conv=1e-3
        )

    # TT+TE+EE+EB
    Ag[6, 2:] = Ag[5, 2:] * Ag[4, 2:] / (Ag[5, 2:] + Ag[4, 2:])
    Ac[6, 2:] = Ac[5, 2:] * Ac[4, 2:] / (Ac[5, 2:] + Ac[4, 2:])

    # EE+EB
    Ag[5, 2:] = Ag[2, 2:] * Ag[4, 2:] / (Ag[2, 2:] + Ag[4, 2:])
    Ac[5, 2:] = Ac[2, 2:] * Ac[4, 2:] / (Ac[2, 2:] + Ac[4, 2:])

    return Ag, Ac


def xtt(est, lmax, rlmin, rlmax, fC, OCT, lfac=''):
    """
    Compute the unnormalized response between lensing potential and amplitude
    modulation from the temperature quadratic estimator.

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
        Theory TT spectrum, with bounds ``0:rlmax``.
    OCT : array_like of float, shape (rlmax + 1,)
        Observed TT spectrum, with bounds ``0:rlmax``.
    lfac : str, optional
        Multiplicative factor. Use ``'k'`` for convergence or ``''`` for
        lensing potential. Default is ``''``.

    Returns
    -------
    Ag : ndarray of float, shape (lmax + 1,)
        Cross normalization, with bounds ``0:lmax``.
    """
    return lib_norm.norm_quad.xtt(est, lmax, rlmin, rlmax, fC, OCT, lfac)


def xeb(est, lmax, rlmin, rlmax, EE, EB, OCE, OCB, BB):
    """
    Compute the response of the reconstructed EB field to another field.

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
    EE : array_like of float, shape (rlmax + 1,)
        Theory EE spectrum, with bounds ``0:rlmax``.
    EB : array_like of float, shape (rlmax + 1,)
        Theory EB spectrum, with bounds ``0:rlmax``.
    OCE : array_like of float, shape (rlmax + 1,)
        Observed EE spectrum, with bounds ``0:rlmax``.
    OCB : array_like of float, shape (rlmax + 1,)
        Observed BB spectrum, with bounds ``0:rlmax``.
    BB : array_like of float, shape (rlmax + 1,)
        Theory BB spectrum, with bounds ``0:rlmax``.

    Returns
    -------
    Aa : ndarray of float, shape (lmax + 1,)
        Polarization-rotation angle normalization, with bounds ``0:lmax``.
    """
    return lib_norm.norm_quad.xeb(est, lmax, rlmin, rlmax, EE, EB, OCE, OCB, BB)
