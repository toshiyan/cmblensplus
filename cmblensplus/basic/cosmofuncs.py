from . import libbasic


def hubble(z, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0., divc=False):
    """
    Compute the expansion rate.

    The returned value is either the expansion rate ``H`` in units of
    ``km/s/Mpc`` or ``H / c`` in units of ``1/Mpc``.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshifts at which the expansion rate is computed.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.
    divc : bool, optional
        If True, divide ``H`` by the speed of light. Default is False.

    Returns
    -------
    Hz : ndarray of float, shape (zn,)
        Expansion rate ``H(z)``, divided by the speed of light if
        ``divc`` is True.
    """
    zn = len(z)
    return libbasic.cosmofuncs.hubble(z, H0, Om, Ov, w0, wa, zn, divc)


def dhubble_dz(z, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0.):
    """
    Compute the derivative of the expansion rate, ``dH(z) / dz``.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshifts at which ``dH/dz`` is computed.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.

    Returns
    -------
    dHdz : ndarray of float, shape (zn,)
        Derivative of the expansion rate, ``dH(z) / dz``.
    """
    zn = len(z)
    return libbasic.cosmofuncs.dhubble_dz(z, H0, Om, Ov, w0, wa, zn)


def dist2z(rz, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0.):
    """
    Compute redshift as a function of comoving distance.

    Parameters
    ----------
    rz : array_like of float, shape (zn,)
        Comoving distance in units of ``Mpc``.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.

    Returns
    -------
    z : ndarray of float, shape (zn,)
        Redshift.
    """
    zn = len(rz)
    return libbasic.cosmofuncs.dist2z(rz, H0, Om, Ov, w0, wa, zn)


def dist_comoving(z, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0.):
    """
    Compute comoving distance as a function of redshift.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshift.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.

    Returns
    -------
    rz : ndarray of float, shape (zn,)
        Comoving distance in units of ``Mpc``.
    """
    zn = len(z)
    return libbasic.cosmofuncs.dist_comoving(z, H0, Om, Ov, w0, wa, zn)


def dist_luminosity(z, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0.):
    """
    Compute luminosity distance as a function of redshift.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshift.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.

    Returns
    -------
    DLz : ndarray of float, shape (zn,)
        Luminosity distance in units of ``Mpc``.
    """
    zn = len(z)
    return libbasic.cosmofuncs.dist_luminosity(z, H0, Om, Ov, w0, wa, zn)


def growth_factor(z, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0., normed=False):
    """
    Compute the analytic linear growth factor ``D(z)`` as a function of
    redshift.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshift.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.
    normed : bool, optional
        If True, normalize the growth factor so that ``D(z=0) = 1``.
        Otherwise, use the normalization in which ``D(z) = a`` in a pure
        matter universe with ``Omega_m(a) = 1``. Default is False.

    Returns
    -------
    Dz : ndarray of float, shape (zn,)
        Growth factor.
    """
    zn = len(z)
    return libbasic.cosmofuncs.growth_factor(
        z, H0, Om, Ov, w0, wa, zn, normed
    )


def growth_rate(z, H0=70., Om=0.3, Ov=0.7, w0=-1., wa=0.):
    """
    Compute the linear growth rate as a function of redshift.

    The growth rate is defined as ``f(z) = d ln D / d ln a``.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshift.
    H0 : float, optional
        Present-day Hubble parameter in units of ``km/s/Mpc``.
        Default is 70.
    Om : float, optional
        Present-day matter density parameter. Default is 0.3.
    Ov : float, optional
        Present-day dark-energy density parameter. Default is 0.7.
    w0 : float, optional
        Dark-energy equation-of-state parameter. Default is -1.
    wa : float, optional
        Time-dependent dark-energy equation-of-state parameter.
        Default is 0.

    Returns
    -------
    fz : ndarray of float, shape (zn,)
        Linear growth rate.
    """
    zn = len(z)
    return libbasic.cosmofuncs.growth_rate(z, H0, Om, Ov, w0, wa, zn)


def nz_gw(z, Cz, Hz, ntype='CH06', dotn0=1e-6, Tobs=3.):
    """
    Compute the redshift distribution of neutron-star binary merger events.

    This returns the event distribution per redshift, ``dN/dz``, at ``z``.

    Parameters
    ----------
    z : float
        Redshift.
    Cz : float
        Comoving distance.
    Hz : float
        Expansion rate.
    ntype : str, optional
        Type of merger-rate functional form. Use ``'CH06'`` for the default
        model or ``'none'`` for no redshift evolution. Default is ``'CH06'``.
    dotn0 : float, optional
        Present-day merger rate. Default is ``1e-6``.
    Tobs : float, optional
        Total observation time. Default is 3.

    Returns
    -------
    nz : float
        Event distribution function at ``z``.
    """
    return libbasic.cosmofuncs.nz_gw(z, Cz, Hz, ntype, dotn0, Tobs)


def drate_dz(z, ntype='CH06'):
    """
    Compute the redshift-dependent merger-rate factor.

    Parameters
    ----------
    z : array_like of float, shape (zn,)
        Redshift.
    ntype : str, optional
        Type of merger-rate functional form. Use ``'CH06'`` for the default
        model or ``'none'`` for no redshift evolution. Default is ``'CH06'``.

    Returns
    -------
    drate : ndarray of float, shape (zn,)
        Redshift-dependent merger-rate factor.
    """
    zn = len(z)
    return libbasic.cosmofuncs.drate_dz(z, zn, ntype)
    