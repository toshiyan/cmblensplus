import numpy
from ._core import lib_bispec


def make_quad_gauss(alm, bst=1):
    r"""
    Return a non-Gaussian alm.

    The corresponding non-Gaussian field is defined as

    .. math::

        \delta^{\rm NL}(\hat{n})
        =
        \delta^{\rm L}(\hat{n})
        +
        \left[\delta^{\rm L}(\hat{n})\right]^2,

    where :math:`\delta^{\rm L}(\hat{n})` is a Gaussian field obtained from
    the input alm.

    Parameters
    ----------
    alm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Input harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``.
    bst : int, optional
        Accuracy-control parameter passed to the backend routine. Default is 1.

    Returns
    -------
    qlm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Output harmonic coefficients of the non-Gaussian field, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_bispec.make_quad_gauss(alm, bst=bst)


def bispec_norm(bp, bstype='equi', bst=2, sL=None, nthreads=0):
    """
    Return the normalization of the binned reduced bispectrum for given
    multipole bins.

    Parameters
    ----------
    bp : array_like of float, shape (bn + 1,)
        Bin edges, with bounds ``0:bn``.
    bstype : str, optional
        Bispectrum configuration. Default is ``'equi'`` for equilateral.
    bst : int, optional
        Accuracy-control parameter, defined as ``bst = nside / lmax``.
        Larger values give more accurate results. Default is 2.
    sL : array_like of int, shape (2,), optional
        Fixed bins for the squeezed configuration, ``b[sL, eL, eL]``.
        If not given, the lowest multipole bin is used.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    norm : ndarray of float, shape (bn,)
        Normalization of the binned reduced bispectrum at each bin, with
        bounds ``0:bn-1``.
    """
    return lib_bispec.bispec_norm(
        bp, bstype=bstype, bst=bst, sL=sL, nthreads=nthreads
    )


def bispec_bin(bp, alm, bstype='equi', bst=2, sL=None, nthreads=0):
    """
    Return the unnormalized binned reduced bispectrum for given multipole bins.

    Parameters
    ----------
    bp : array_like of float, shape (bn + 1,)
        Bin edges, with bounds ``0:bn``.
    alm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Input harmonic coefficients, with bounds ``(0:lmax, 0:lmax)``.
    bstype : str, optional
        Bispectrum configuration. Default is ``'equi'`` for equilateral.
    bst : int, optional
        Accuracy-control parameter, defined as ``bst = nside / lmax``.
        Larger values give more accurate results. Default is 2.
    sL : array_like of int, shape (2,), optional
        Fixed bins for the squeezed configuration, ``b[sL, eL, eL]``.
        If not given, the lowest multipole bin is used.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    bis : ndarray of float, shape (bn,)
        Unnormalized binned reduced bispectrum at each bin, with bounds
        ``0:bn-1``.
    """
    return lib_bispec.bispec_bin(
        bp, alm, bstype, bst, sL, nthreads=nthreads
    )


def xbispec_bin(bp, alm, bstype='equi', bst=2, sL=None, nthreads=0):
    """
    Return the unnormalized binned reduced cross-bispectrum for given
    multipole bins.

    Parameters
    ----------
    bp : array_like of float, shape (bn + 1,)
        Bin edges, with bounds ``0:bn``.
    alm : ndarray of complex, shape (n, lmax + 1, lmax + 1)
        Input harmonic coefficients, with bounds
        ``(0:n-1, 0:lmax, 0:lmax)``.
    bstype : str, optional
        Bispectrum configuration. Default is ``'equi'`` for equilateral.
    bst : int, optional
        Accuracy-control parameter, defined as ``bst = nside / lmax``.
        Larger values give more accurate results. Default is 2.
    sL : array_like of int, shape (2,), optional
        Fixed bins for the squeezed configuration, ``b[sL, eL, eL]``.
        If not given, the lowest multipole bin is used.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    bis : ndarray of float, shape (bn,)
        Unnormalized binned reduced cross-bispectrum at each bin, with bounds
        ``0:bn-1``.
    """
    return lib_bispec.xbispec_bin(
        bp, alm, bstype, bst, sL, nthreads=nthreads
    )