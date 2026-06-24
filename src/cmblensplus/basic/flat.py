from . import libbasic


def alxy(qest, qtype, lmax, rlmin, rlmax, fC, W1, W2, gln=100, gle=1e-14, lxcut=0):
    r"""
    Compute the flat-sky quadratic-estimator normalization.

    Notes
    -----
    This routine interpolates the input :math:`C_\ell` at non-integer multipoles by using ``Cl[int(ell)]``. 
    This may lead to a small discrepancy between the normalization computed by this routine and the FFT-based
    normalization, which uses linear interpolation.

    The FFT-based normalization is recommended when normalizing simulation results.

    Parameters
    ----------
    qest : {'TT', 'TE', 'TB', 'EE', 'EB', 'BB'}
        Estimator combination.
    qtype : str
        Estimator type, for example ``'lensing'`` or ``'patchytau'``.
    lmax : int
        Maximum multipole of the output normalization.
    rlmin : int
        Minimum input CMB multipole used for reconstruction.
    rlmax : int
        Maximum input CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        Power spectrum in the numerator, with bounds ``0:rlmax``.
    W1 : array_like of float, shape (rlmax + 1,)
        Inverse observed power spectrum for the first field, with bounds ``0:rlmax``.
    W2 : array_like of float, shape (rlmax + 1,)
        Inverse observed power spectrum for the second field, with bounds ``0:rlmax``.
    gln : int, optional
        Number of Gauss-Legendre integration points. Default is 100.
    gle : float, optional
        Convergence parameter for the Gauss-Legendre integration. Default is ``1e-14``.
    lxcut : int, optional
        Multipole cut in the x direction, :math:`|\ell_x| <` ``lxcut``. Default is 0.

    Returns
    -------
    Ag : ndarray of float, shape (lmax + 1,)
        Normalization for the even estimator pair.
    Ac : ndarray of float, shape (lmax + 1,)
        Normalization for the odd estimator pair.
    """
    return libbasic.flat.alxy(
        qest, qtype, lmax, rlmin, rlmax, fC, W1, W2, gln, gle, lxcut
    )


def alxy_asym(qest, qtype, lmax, rlmin, rlmax, fC, AA, BB, AB, gln=100, gle=1e-14, lxcut=0):
    r"""
    Compute the asymmetric flat-sky quadratic-estimator normalization.

    Parameters
    ----------
    qest : {'TT', 'TE', 'TB', 'EE', 'EB', 'BB'}
        Estimator combination.
    qtype : str
        Estimator type, for example ``'lensing'`` or ``'patchytau'``.
    lmax : int
        Maximum multipole of the output normalization.
    rlmin : int
        Minimum input CMB multipole used for reconstruction.
    rlmax : int
        Maximum input CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        Power spectrum in the numerator, with bounds ``0:rlmax``.
    AA : array_like of float, shape (rlmax + 1,)
        Observed auto spectrum of the first field, with bounds ``0:rlmax``.
    BB : array_like of float, shape (rlmax + 1,)
        Observed auto spectrum of the second field, with bounds ``0:rlmax``.
    AB : array_like of float, shape (rlmax + 1,)
        Observed cross spectrum between the first and second fields, with bounds ``0:rlmax``.
    gln : int, optional
        Number of Gauss-Legendre integration points. Default is 100.
    gle : float, optional
        Convergence parameter for the Gauss-Legendre integration. Default is ``1e-14``.
    lxcut : int, optional
        Multipole cut in the x direction, :math:`|\ell_x| <` ``lxcut``. Default is 0.

    Returns
    -------
    Ag : ndarray of float, shape (lmax + 1,)
        Normalization for the even estimator pair.
    Ac : ndarray of float, shape (lmax + 1,)
        Normalization for the odd estimator pair.
    """
    return libbasic.flat.alxy_asym(
        qest, qtype, lmax, rlmin, rlmax, fC, AA, BB, AB, gln, gle, lxcut
    )


def bbxy(lmax, rlmin, rlmax, XX, YY, weight='lensing', gln=100, gle=1e-14):
    """
    Compute the flat-sky B-mode power spectrum induced by a modulation field to E-mode or T.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output normalization.
    rlmin : int
        Minimum input multipole used in the calculation.
    rlmax : int
        Maximum input multipole used in the calculation.
    XX : array_like of float, shape (rlmax + 1,)
        Modulation field power spectrum such as lensing potential power spectrum, with bounds ``0:rlmax``.
    YY : array_like of float, shape (rlmax + 1,)
        CMB E or T auto power spectrum, with bounds ``0:rlmax``.
    weight : str, optional
        Weight type. Default is ``'lensing'``.
    gln : int, optional
        Number of Gauss-Legendre integration points. Default is 100.
    gle : float, optional
        Convergence parameter for the Gauss-Legendre integration. Default is ``1e-14``.

    Returns
    -------
    BB : ndarray of float, shape (lmax + 1,)
        Normalization, with bounds ``0:lmax``.
    """
    return libbasic.flat.bbxy(lmax, rlmin, rlmax, XX, YY, weight, gln, gle)
