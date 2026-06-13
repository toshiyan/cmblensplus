from . import libflatsky


def bispec_norm(nx, ny, D, bp, dbin_max=-1, bn=1):
    """
    Return the normalization of the binned bispectrum estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    bp : array_like of float, shape (bn + 1,)
        Multipole bin edges.
    dbin_max : int, optional
        Maximum bin separation used in the bispectrum calculation. If -1,
        it is set to ``bn``. Default is -1.
    bn : int, optional
        Number of multipole bins. This value is overwritten by
        ``len(bp) - 1``. Default is 1.

    Returns
    -------
    norm : ndarray of float
        Normalization of the binned bispectrum estimator.
    """
    bn = len(bp) - 1
    if dbin_max == -1:
        dbin_max = bn
    return libflatsky.bispec.bispec_norm(nx, ny, D, bp, dbin_max, bn)


def bispec_bin(kmap, bp, kn=1, bn=1, nx=0, ny=0, dbin_max=-1):
    """
    Return the binned bispectrum estimator.

    Parameters
    ----------
    kmap : ndarray of complex, shape (kn, bn, nx, ny)
        Fourier-space maps used for the bispectrum calculation.
    bp : array_like of float, shape (bn + 1,)
        Multipole bin edges.
    kn : int, optional
        Number of input map sets. This value is overwritten from ``kmap``.
        Default is 1.
    bn : int, optional
        Number of multipole bins. This value is overwritten by
        ``len(bp) - 1``. Default is 1.
    nx : int, optional
        Number of Fourier grid points along the x direction. This value is
        overwritten from ``kmap``. Default is 0.
    ny : int, optional
        Number of Fourier grid points along the y direction. This value is
        overwritten from ``kmap``. Default is 0.
    dbin_max : int, optional
        Maximum bin separation used in the bispectrum calculation. If -1,
        it is set to ``bn``. Default is -1.

    Returns
    -------
    bispec : ndarray of float
        Binned bispectrum estimator.
    """
    bn = len(bp) - 1
    kn = len(kmap[:, 0, 0, 0])
    nx = len(kmap[0, 0, :, 0])
    ny = len(kmap[0, 0, 0, :])
    if dbin_max == -1:
        dbin_max = bn
    return libflatsky.bispec.bispec_bin(kn, bn, nx, ny, kmap, bp, dbin_max)


def binfilter(nx, ny, D, bp, bn=1):
    """
    Return binary masks for multipole bins on two-dimensional Fourier grids.

    The mask is 1 inside the corresponding multipole bin and 0 otherwise.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    bp : array_like of float, shape (bn + 1,)
        Multipole bin edges.
    bn : int, optional
        Number of multipole bins. This value is overwritten by
        ``len(bp) - 1``. Default is 1.

    Returns
    -------
    bf : ndarray of float, shape (bn, nx, ny)
        Binary mask for each multipole bin.
    """
    bn = len(bp) - 1
    return libflatsky.bispec.binfilter(nx, ny, D, bp, bn)


def bispec_norm_1d(nx, ny, D, bfs, bn=1):
    """
    Return the normalization of the one-dimensional binned bispectrum
    estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    bfs : ndarray of float, shape (3, bn, nx, ny)
        Multipole-bin masks on two-dimensional grids, usually obtained from
        :func:`binfilter`.
    bn : int, optional
        Number of multipole bins. This value is overwritten from ``bfs``.
        Default is 1.

    Returns
    -------
    bnorm : ndarray of float, shape (bn,)
        Normalization of the one-dimensional binned bispectrum estimator at
        each multipole bin.
    """
    bn = len(bfs[0, :, 0, 0])
    return libflatsky.bispec.bispec_norm_1d(nx, ny, D, bfs, bn)


def bispec_bin_1d(nx, ny, D, bfs, bnorm, alm, bn=1):
    """
    Return the one-dimensional binned bispectrum estimator.

    Parameters
    ----------
    nx : int
        Number of Fourier grid points along the x direction.
    ny : int
        Number of Fourier grid points along the y direction.
    D : array_like of float, shape (2,)
        Map side lengths, equivalent to ``dLx / (2*pi)`` and
        ``dLy / (2*pi)``.
    bfs : ndarray of float, shape (3, bn, nx, ny)
        Multipole-bin masks on two-dimensional grids, usually obtained from
        :func:`binfilter`.
    bnorm : array_like of float, shape (bn,)
        Normalization of the one-dimensional binned bispectrum estimator at
        each multipole bin.
    alm : ndarray of complex, shape (3, nx, ny)
        Fourier modes for each leg of the bispectrum.
    bn : int, optional
        Number of multipole bins. This value is overwritten from ``bfs``.
        Default is 1.

    Returns
    -------
    bispec : ndarray of float, shape (bn,)
        One-dimensional binned bispectrum at each multipole bin.
    """
    bn = len(bfs[0, :, 0, 0])
    return libflatsky.bispec.bispec_bin_1d(nx, ny, D, bfs, bnorm, alm, bn)
    