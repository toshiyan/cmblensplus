from . import libbasic
import numpy


def cl_flat(cpmodel, z, dz, zs, lmax, k, pk0, pktype='T12', cltype='kk', dNdz=None, wdel=None):
    r"""
    Compute the flat-sky lensing power spectrum analytically.

    Parameters
    ----------
    cpmodel : str
        Cosmological parameter model. Supported values include ``'model0'``,
        ``'modelw'``, and ``'modelp'``.
    z : array_like of float, shape (zn,)
        Redshift points for the redshift integral.
    dz : array_like of float, shape (zn,)
        Redshift intervals.
    zs : array_like of float, shape (2,)
        Source redshifts.
    lmax : int
        Maximum multipole of the output power spectrum.
    k : array_like of float, shape (kn,)
        Wavenumbers for the matter power spectrum in units of
        :math:`h/{\rm Mpc}`.
    pk0 : array_like of float, shape (kn,)
        Linear matter power spectrum at ``z = 0`` in units of
        :math:`{\rm Mpc}^3/h^3`.
    pktype : str, optional
        Fitting formula for the matter power spectrum. Supported values include
        ``'Lin'``, ``'S02'``, and ``'T12'``. Default is ``'T12'``.
    cltype : str, optional
        Power-spectrum type. Supported values include ``'kk'``, ``'gk'``, and
        ``'gg'``. Default is ``'kk'``.
    dNdz : array_like of float, shape (zn,), optional
        Redshift distribution of galaxies. Used only when ``cltype`` includes
        ``'g'``.
    wdel : array_like of float, shape (zn, lmax + 1), optional
        Modified chi-kernel function for z-cleaning from ``l = 0`` to
        ``lmax``.

    Returns
    -------
    cl : ndarray of float, shape (lmax + 1,)
        Power spectrum from LSS contributions.
    """
    zn, kn = len(z), len(k)
    if dNdz is None:
        dNdz = z * 0.
    if wdel is None:
        wdel = numpy.zeros((zn, lmax + 1))
    if len(dNdz) != zn:
        raise ValueError('len(dNdz) should match len(z)')
    return libbasic.bispec.cl_flat(cpmodel, z, dz, zs, lmax, k, pk0, zn, kn, pktype, cltype, dNdz, wdel)


def bispeclens(shap, cpmodel, model, z, dz, zs, lmin, lmax, k, pk0, lan=0., kan=0., pktype='T12', ltype='', btype='kkk', dNdz=None, wdel=None):
    r"""
    Compute the lensing bispectrum analytically.

    Parameters
    ----------
    shap : str
        Bispectrum shape. Supported values include ``'equi'``, ``'fold'``,
        ``'sque'``, and ``'isos'``.
    cpmodel : str
        Cosmological parameter model. Supported values include ``'model0'``,
        ``'modelw'``, and ``'modelp'``.
    model : str
        Matter-bispectrum fitting formula. Supported values include
        ``'LN'``, ``'SC'``, ``'GM'``, ``'3B'``, and ``'RT'``.
    z : array_like of float, shape (zn,)
        Redshift points for the redshift integral.
    dz : array_like of float, shape (zn,)
        Redshift intervals.
    zs : array_like of float, shape (3,)
        Source redshifts.
    lmin : int
        Minimum multipole of the bispectrum.
    lmax : int
        Maximum multipole of the bispectrum.
    k : array_like of float, shape (kn,)
        Wavenumbers for the matter power spectrum in units of
        :math:`h/{\rm Mpc}`.
    pk0 : array_like of float, shape (kn,)
        Linear matter power spectrum at ``z = 0`` in units of
        :math:`{\rm Mpc}^3/h^3`.
    lan : float, optional
        Modified-gravity extension parameter. Default is 0.
    kan : float, optional
        Modified-gravity extension parameter. Default is 0.
    pktype : str, optional
        Fitting formula for the matter power spectrum. Supported values include
        ``'Lin'``, ``'S02'``, and ``'T12'``. Default is ``'T12'``.
    ltype : str, optional
        Full-sky correction option. Use ``'full'`` to include the correction.
        Default is ``''``.
    btype : str, optional
        Bispectrum type. Supported values include ``'kkk'``, ``'gkk'``,
        ``'ggk'``, and ``'ggg'``. Default is ``'kkk'``.
    dNdz : array_like of float, shape (zn,), optional
        Redshift distribution of galaxies. Used only when ``btype`` includes
        ``'g'``.
    wdel : array_like of float, shape (zn, lmax + 1), optional
        Modified chi-kernel function for z-cleaning from ``l = 0`` to
        ``lmax``.

    Returns
    -------
    bl0 : ndarray of float, shape (lmax + 1,)
        Lensing bispectrum from LSS contributions.
    bl1 : ndarray of float, shape (lmax + 1,)
        Lensing bispectrum from post-Born contributions.
    """
    zn, kn = len(z), len(k)
    if dNdz is None:
        dNdz = z * 0.
    if wdel is None:
        wdel = numpy.zeros((zn, lmax + 1))
    if len(dNdz) != zn:
        raise SystemExit('ERROR in "bispeclens": size of dNdz should be zn')
    return libbasic.bispec.bispeclens(shap, cpmodel, model, z, dz, zs, lmin, lmax, k, pk0, lan, kan, zn, kn, pktype, ltype, btype, dNdz, wdel)


def bispeclens_bin(shap, cpmodel, model, z, dz, zs, lmin, lmax, bn, k, pk0, lan=0., kan=0., pktype='T12', btype='kkk', dNdz=None, wdel=None):
    r"""
    Compute the binned lensing bispectrum analytically.

    Parameters
    ----------
    shap : str
        Bispectrum shape. Supported values include ``'equi'``, ``'fold'``,
        ``'sque'``, and ``'isos'``.
    cpmodel : str
        Cosmological parameter model. Supported values include ``'model0'``,
        ``'modelw'``, and ``'modelp'``.
    model : str
        Matter-bispectrum fitting formula. Supported values include
        ``'LN'``, ``'SC'``, ``'GM'``, ``'3B'``, and ``'RT'``.
    z : array_like of float, shape (zn,)
        Redshift points for the redshift integral.
    dz : array_like of float, shape (zn,)
        Redshift intervals.
    zs : array_like of float, shape (3,)
        Source redshifts.
    lmin : int
        Minimum multipole of the bispectrum.
    lmax : int
        Maximum multipole of the bispectrum.
    bn : int
        Number of multipole bins.
    k : array_like of float, shape (kn,)
        Wavenumbers for the matter power spectrum in units of
        :math:`h/{\rm Mpc}`.
    pk0 : array_like of float, shape (kn,)
        Linear matter power spectrum at ``z = 0`` in units of
        :math:`{\rm Mpc}^3/h^3`.
    lan : float, optional
        Modified-gravity extension parameter. Default is 0.
    kan : float, optional
        Modified-gravity extension parameter. Default is 0.
    pktype : str, optional
        Fitting formula for the matter power spectrum. Supported values include
        ``'Lin'``, ``'S02'``, and ``'T12'``. Default is ``'T12'``.
    btype : str, optional
        Bispectrum type. Supported values include ``'kkk'``, ``'gkk'``,
        ``'ggk'``, and ``'ggg'``. Default is ``'kkk'``.
    dNdz : array_like of float, shape (zn,), optional
        Redshift distribution of galaxies. Used only when ``btype`` includes
        ``'g'``.
    wdel : array_like of float, shape (zn, lmax + 1), optional
        Modified chi-kernel function for z-cleaning from ``l = 0`` to
        ``lmax``.

    Returns
    -------
    bc : ndarray of float, shape (bn,)
        Multipole bin centers.
    bl0 : ndarray of float, shape (bn,)
        Binned lensing bispectrum from LSS contributions.
    bl1 : ndarray of float, shape (bn,)
        Binned lensing bispectrum from post-Born contributions.
    """
    zn, kn = len(z), len(k)
    if dNdz is None:
        dNdz = z * 0.
    if wdel is None:
        wdel = numpy.zeros((zn, lmax + 1))
    if len(dNdz) != zn:
        raise SystemExit('ERROR in "bispeclens_bin": size of dNdz should be zn')
    return libbasic.bispec.bispeclens_bin(shap, cpmodel, model, z, dz, zn, zs, lmin, lmax, bn, k, pk0, kn, lan, kan, pktype, btype, dNdz, wdel)


def bispeclens_snr(cpmodel, model, z, dz, zs, lmin, lmax, cl, k, pk0, pktype='T12', btype='kkk', dNdz=None, cgg=None, ro=100, wdel=None):
    r"""
    Compute the signal-to-noise ratio of the lensing bispectrum analytically.

    Parameters
    ----------
    cpmodel : str
        Cosmological parameter model. Supported values include ``'model0'``,
        ``'modelw'``, and ``'modelp'``.
    model : str
        Matter-bispectrum fitting formula. Supported values include
        ``'LN'``, ``'SC'``, ``'GM'``, ``'3B'``, and ``'RT'``.
    z : array_like of float, shape (zn,)
        Redshift points for the redshift integral.
    dz : array_like of float, shape (zn,)
        Redshift intervals.
    zs : array_like of float, shape (3,)
        Source redshifts.
    lmin : int
        Minimum multipole of the bispectrum.
    lmax : int
        Maximum multipole of the bispectrum.
    cl : array_like of float, shape (lmax + 1,)
        Observed angular power spectrum from ``l = 0`` to ``lmax``.
    k : array_like of float, shape (kn,)
        Wavenumbers for the matter power spectrum in units of
        :math:`h/{\rm Mpc}`.
    pk0 : array_like of float, shape (kn,)
        Linear matter power spectrum at ``z = 0`` in units of
        :math:`{\rm Mpc}^3/h^3`.
    pktype : str, optional
        Fitting formula for the matter power spectrum. Supported values include
        ``'Lin'``, ``'S02'``, and ``'T12'``. Default is ``'T12'``.
    btype : str, optional
        Bispectrum type. Supported values include ``'kkk'``, ``'gkk'``,
        ``'ggk'``, and ``'ggg'``. Default is ``'kkk'``.
    dNdz : array_like of float, shape (zn,), optional
        Redshift distribution of galaxies. Used only when ``btype`` includes
        ``'g'``.
    cgg : array_like of float, shape (lmax + 1,), optional
        Observed galaxy spectrum. If not given, a zero spectrum is used.
    ro : int, optional
        Progress-output interval in multipoles. Default is 100.
    wdel : array_like of float, shape (zn, lmax + 1), optional
        Modified chi-kernel function for z-cleaning from ``l = 0`` to
        ``lmax``.

    Returns
    -------
    snr : ndarray of float, shape (2,)
        Total signal-to-noise ratio and LSS-only signal-to-noise ratio.
    """
    zn, kn = len(z), len(k)
    if dNdz is None:
        dNdz = z * 0.
    if len(dNdz) != zn:
        raise SystemExit('ERROR in "bispeclens_snr": size of dNdz should be zn')
    if wdel is None:
        wdel = numpy.zeros((zn, lmax + 1))
    if cgg is None:
        cgg = cl * 0.
    if len(cgg) != lmax + 1:
        raise SystemExit('ERROR in "bispeclens_snr": size of cgg should be lmax+1')
    if lmin < 1:
        raise SystemExit('ERROR in "bispeclens_snr": lmin should be >=1')
    return libbasic.bispec.bispeclens_snr(cpmodel, model, z, dz, zn, zs, lmin, lmax, cl, k, pk0, kn, pktype, btype, dNdz, cgg, ro, wdel)


def bispeclens_gauss_bin(shap, bn, lmin, lmax, cl):
    """
    Compute the binned bispectrum analytically for the quadratic Gaussian
    model.

    Parameters
    ----------
    shap : str
        Bispectrum shape. Supported values include ``'equi'``, ``'fold'``,
        ``'sque'``, and ``'isos'``.
    bn : int
        Number of multipole bins.
    lmin : int
        Minimum multipole of the bispectrum.
    lmax : int
        Maximum multipole of the bispectrum.
    cl : array_like of float, shape (lmax + 1,)
        Power spectrum from ``l = 0`` to ``lmax``.

    Returns
    -------
    bc : ndarray of float, shape (bn,)
        Multipole bin centers.
    bl : ndarray of float, shape (bn,)
        Binned bispectrum.
    """
    if lmin < 1:
        raise SystemExit('ERROR in "bispeclens_gauss_bin": lmin should be >=1')
    return libbasic.bispec.bispeclens_gauss_bin(shap, bn, lmin, lmax, cl)


def zpoints(zmin, zmax, zn, zspace=1):
    """
    Precompute interpolation points for redshift.

    Parameters
    ----------
    zmin : float
        Minimum redshift.
    zmax : float
        Maximum redshift.
    zn : int
        Number of redshift points.
    zspace : int, optional
        Type of spacing. Use 0 for linear spacing and 1 for Gauss-Legendre
        spacing. Default is 1.

    Returns
    -------
    z : ndarray of float, shape (zn,)
        Redshift points.
    dz : ndarray of float, shape (zn,)
        Redshift intervals.
    """
    if not zspace in [0, 1]:
        raise SystemExit('ERROR in "zpoints": zspace should be 0 or 1')
    return libbasic.bispec.zpoints(zmin, zmax, zn, zspace)


def skewspeclens(
    cpmodel, model, z, dz, zs, ols, lmin, lmax, k, pk0, theta=0.0,
    pktype='T12', btype='kkk', pb=True, Om=0.3, H0=70., w0=-1., wa=0.,
    mnu=0.06, ns=0.965, verbose=True, dNdz=None, wdel=None
):
    r"""
    Compute the skew spectrum using a matter-bispectrum fitting formula.

    Parameters
    ----------
    cpmodel : str
        Cosmological parameter model. Supported values include ``'model0'``,
        ``'modelw'``, ``'modelp'``, and ``'input'``.
    model : str
        Matter-bispectrum fitting formula. Supported values include
        ``'LN'``, ``'SC'``, ``'GM'``, ``'3B'``, and ``'RT'``.
    z : array_like of float, shape (zn,)
        Redshift points for the redshift integral.
    dz : array_like of float, shape (zn,)
        Redshift intervals.
    zs : array_like of float, shape (2,)
        Source redshifts. The second source redshift is used for the squared
        map.
    ols : array_like of int, shape (bn,)
        Output multipoles to be computed for the skew spectrum.
    lmin : int
        Minimum multipole of alms included in the skew spectrum.
    lmax : int
        Maximum multipole of alms included in the skew spectrum.
    k : array_like of float, shape (kn,)
        Wavenumbers for the matter power spectrum in units of
        :math:`h/{\rm Mpc}`.
    pk0 : array_like of float, shape (kn,)
        Linear matter power spectrum at ``z = 0`` in units of
        :math:`{\rm Mpc}^3/h^3`.
    theta : float, optional
        Kappa-map resolution in arcmin. Default is 0.0.
    pktype : str, optional
        Fitting formula for the matter power spectrum. Supported values include
        ``'Lin'``, ``'S02'``, and ``'T12'``. Default is ``'T12'``.
    btype : str, optional
        Bispectrum type. Supported values include ``'kkk'``, ``'gkk'``,
        ``'kgg'``, and ``'ggg'``. Default is ``'kkk'``.
    pb : bool, optional
        Whether to include the post-Born correction. Default is True.
    Om : float, optional
        Matter density parameter used when ``cpmodel='input'``. Default is
        0.3.
    H0 : float, optional
        Hubble constant used when ``cpmodel='input'``. Default is 70.
    w0 : float, optional
        Dark-energy equation-of-state parameter used when ``cpmodel='input'``.
        Default is -1.
    wa : float, optional
        Dark-energy equation-of-state parameter used when ``cpmodel='input'``.
        Default is 0.
    mnu : float, optional
        Sum of neutrino masses used when ``cpmodel='input'``. Default is 0.06.
    ns : float, optional
        Scalar spectral index used when ``cpmodel='input'``. Default is 0.965.
    verbose : bool, optional
        Whether to output messages. Default is True.
    dNdz : array_like of float, shape (zn,), optional
        Redshift distribution of galaxies. Used only when ``btype`` includes
        ``'g'``.
    wdel : array_like of float, shape (zn, lmax + 1), optional
        Modified chi-kernel function for z-cleaning from ``l = 0`` to
        ``lmax``.

    Returns
    -------
    skew : ndarray of float, shape (3, 2, bn)
        Skew spectra ``S0``, ``S1``, and ``S2`` from LSS and post-Born
        contributions, separately.
    """
    bn, zn, kn = len(ols), len(z), len(k)
    if dNdz is None:
        dNdz = z * 0.
    if wdel is None:
        wdel = numpy.zeros((zn, lmax + 1))
    if len(dNdz) != zn:
        raise SystemExit('ERROR in "skewspeclens": len(dNdz) should be len(z)')
    return libbasic.bispec.skewspeclens(
        cpmodel, model, z, dz, zs, ols, lmin, lmax, k, pk0, bn, zn, kn, theta,
        pktype, btype, pb, Om, H0, w0, wa, mnu, ns, verbose, dNdz, wdel
    )
