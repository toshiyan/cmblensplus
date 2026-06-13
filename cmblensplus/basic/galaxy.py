from . import libbasic


def dndz_sf(z, a, b, zm=0, z0=0):
    """
    Return a model galaxy redshift distribution.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshifts at which ``dndz`` is returned.
    a : float
        Shape parameter of the Schechter-like galaxy distribution.
    b : float
        Shape parameter of the Schechter-like galaxy distribution.
    zm : float, optional
        Mean redshift. Default is 0.
    z0 : float, optional
        Parameter related to ``zm``. Default is 0. Either ``zm`` or ``z0``
        has to be specified.

    Returns
    -------
    dndz : ndarray of float, shape (zn,)
        Galaxy redshift distribution.
    """
    zn = len(z)
    if zm <= 0 and z0 <= 0:
        raise SystemExit('ERROR in "dndz_sf": both zm and z0 <=0')
    return libbasic.galaxy.dndz_sf(zn, z, a, b, zm, z0)


def photoz_error(z, zi, zn=None, sigma=0.03, zbias=0.):
    """
    Return the photo-z error function for a redshift distribution.

    The returned function is multiplied by the original galaxy redshift
    distribution. See Eq. (13) of arXiv:1103.1118 for details.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshifts at which the photo-z error function is returned.
    zi : array_like of float, shape (2,)
        Redshift-bin edges.
    zn : int, optional
        Number of redshift samples. If not given, ``len(z)`` is used.
    sigma : float, optional
        Photo-z error parameter, used as ``sigma * (1 + z)``.
        Default is 0.03.
    zbias : float, optional
        Photo-z mean bias. Default is 0.

    Returns
    -------
    pz : ndarray of float, shape (zn,)
        Photo-z error function.
    """
    if zn is None:
        zn = len(z)
    return libbasic.galaxy.photoz_error(zn, z, zi, sigma, zbias)


def zbin(zn, a, b, zm=0, z0=0, verbose=False):
    """
    Compute redshift-bin intervals with equal galaxy counts.

    Parameters
    ----------
    zn : int
        Number of redshift bins.
    a : float
        Shape parameter of the Schechter-like galaxy distribution.
    b : float
        Shape parameter of the Schechter-like galaxy distribution.
    zm : float, optional
        Mean redshift. Default is 0.
    z0 : float, optional
        Parameter related to ``zm``. Default is 0. Either ``zm`` or ``z0``
        has to be specified.
    verbose : bool, optional
        Whether to output messages. Default is False.

    Returns
    -------
    zb : ndarray of float, shape (zn + 1,)
        Redshift-bin intervals.
    """
    return libbasic.galaxy.zbin(zn, a, b, zm, z0, verbose)


def frac(zn, zb, a, b, zm, zbias=0.0, sigma=0.0, verbose=False):
    r"""
    Compute the fraction of galaxies in each redshift bin.

    Parameters
    ----------
    zn : int
        Number of redshift bins.
    zb : array_like of float, shape (zn + 1,)
        Redshift-bin intervals.
    a : float
        Shape parameter of the Schechter-like galaxy distribution.
    b : float
        Shape parameter of the Schechter-like galaxy distribution.
    zm : float
        Mean redshift.
    zbias : float, optional
        Constant bias to the true photo-z. Default is 0.
    sigma : float, optional
        Photo-z uncertainty. Default is 0.
    verbose : bool, optional
        Whether to output messages. Default is False.

    Returns
    -------
    nfrac : ndarray of float, shape (zn,)
        Fraction of galaxies in each bin, defined by

        .. math::

            \frac{\int_{z_i}^{z_{i+1}} dz\,N(z)}
                 {\int dz\,N(z)}.
    """
    return libbasic.galaxy.frac(zn, zb, a, b, zm, zbias, sigma, verbose)
