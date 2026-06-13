from . import libflatsky
import numpy as np


def dft1d(map0, nx, ny, npix, D, trans, map1):
    """
    Perform a discrete Fourier transform for a one-dimensional array storing
    data on a two-dimensional grid.

    Parameters
    ----------
    map0 : ndarray of complex, shape (npix,)
        Input data on a two-dimensional grid, with bounds ``0:npix-1``.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    npix : int
        Total number of grid points, ``npix = nx * ny``.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    trans : int
        Transform direction. Use 1 for map to Fourier transform and -1 for
        Fourier to map transform.
    map1 : ndarray of complex, shape (npix,)
        Output array for the transformed data, with bounds ``0:npix-1``.

    Returns
    -------
    map1 : ndarray of complex, shape (npix,)
        Transformed data on the two-dimensional grid, with bounds
        ``0:npix-1``.

    Examples
    --------
    >>> map1 = flatsky.ffttools.dft1d(map0, nx, ny, npix, D, trans, map1)
    """
    return libflatsky.ffttools.dft1d(map0, nx, ny, npix, D, trans, map1)


def dft2d(map0, nx, ny, D, trans):
    """
    Perform a discrete Fourier transform for a two-dimensional complex array.

    Parameters
    ----------
    map0 : ndarray of complex, shape (nx, ny)
        Input data on a two-dimensional grid, with bounds ``(0:nx-1, 0:ny-1)``.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    trans : int
        Transform direction. Use 1 for map to Fourier transform and -1 for
        Fourier to map transform.

    Returns
    -------
    map1 : ndarray of complex, shape (nx, ny)
        Transformed data on the two-dimensional grid, with bounds
        ``(0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> map1 = flatsky.ffttools.dft2d(map0, nx, ny, D, trans)
    """
    return libflatsky.ffttools.dft2d(map0, nx, ny, D, trans)


def dft2dr(map0, nx, ny, D, trans):
    """
    Perform a discrete Fourier transform for a two-dimensional real array.

    Parameters
    ----------
    map0 : ndarray of float, shape (nx, ny)
        Input data on a two-dimensional grid, with bounds ``(0:nx-1, 0:ny-1)``.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    trans : int
        Transform direction. Use 1 for map to Fourier transform and -1 for
        Fourier to map transform.

    Returns
    -------
    map1 : ndarray of float, shape (nx, ny)
        Transformed data on the two-dimensional grid, with bounds
        ``(0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> map1 = flatsky.ffttools.dft2dr(map0, nx, ny, D, trans)
    """
    return libflatsky.ffttools.dft2dr(map0, nx, ny, D, trans)


def dft2drc(map0, nx, ny, D, trans):
    """
    Perform a real-to-complex discrete Fourier transform for a two-dimensional
    array.

    Parameters
    ----------
    map0 : ndarray of float, shape (nx, ny)
        Input real data on a two-dimensional grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    trans : int
        Transform direction. Use 1 for map to Fourier transform and -1 for
        Fourier to map transform.

    Returns
    -------
    map1 : ndarray of complex, shape (nx, ny)
        Transformed data on the two-dimensional grid, with bounds
        ``(0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> map1 = flatsky.ffttools.dft2drc(map0, nx, ny, D, trans)
    """
    return libflatsky.ffttools.dft2drc(map0, nx, ny, D, trans)


def dft2dcr(map0, nx, ny, D, trans):
    """
    Perform a complex-to-real discrete Fourier transform for a two-dimensional
    array.

    Parameters
    ----------
    map0 : ndarray of complex, shape (nx, ny)
        Input complex data on a two-dimensional grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    trans : int
        Transform direction. Use 1 for map to Fourier transform and -1 for
        Fourier to map transform.

    Returns
    -------
    map1 : ndarray of float, shape (nx, ny)
        Transformed data on the two-dimensional grid, with bounds
        ``(0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> map1 = flatsky.ffttools.dft2dcr(map0, nx, ny, D, trans)
    """
    return libflatsky.ffttools.dft2dcr(map0, nx, ny, D, trans)


def dft2dpol(nx, ny, D, Q, U):
    """
    Perform a spin-2 discrete Fourier transform for polarization maps on a
    two-dimensional grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    Q : ndarray of float, shape (nx, ny)
        Q map on the two-dimensional grid, with bounds ``(0:nx-1, 0:ny-1)``.
    U : ndarray of float, shape (nx, ny)
        U map on the two-dimensional grid, with bounds ``(0:nx-1, 0:ny-1)``.

    Returns
    -------
    E : ndarray of complex, shape (nx, ny)
        E-mode map on the two-dimensional Fourier grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    B : ndarray of complex, shape (nx, ny)
        B-mode map on the two-dimensional Fourier grid, with bounds
        ``(0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> E, B = flatsky.ffttools.dft2dpol(nx, ny, D, Q, U)
    """
    return libflatsky.ffttools.dft2dpol(nx, ny, D, Q, U)


def idft2dpol(nx, ny, D, E, B):
    """
    Perform a spin-2 inverse discrete Fourier transform for polarization maps
    on a two-dimensional grid.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    E : ndarray of complex, shape (nx, ny)
        E-mode map on the two-dimensional Fourier grid, with bounds
        ``(0:nx-1, 0:ny-1)``.
    B : ndarray of complex, shape (nx, ny)
        B-mode map on the two-dimensional Fourier grid, with bounds
        ``(0:nx-1, 0:ny-1)``.

    Returns
    -------
    Q : ndarray of float, shape (nx, ny)
        Q map on the two-dimensional grid, with bounds ``(0:nx-1, 0:ny-1)``.
    U : ndarray of float, shape (nx, ny)
        U map on the two-dimensional grid, with bounds ``(0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> Q, U = flatsky.ffttools.idft2dpol(nx, ny, D, E, B)
    """
    return libflatsky.ffttools.idft2dpol(nx, ny, D, E, B)


def eb_separate(nx, ny, D, QU, W, Wd=None):
    """
    Compute Smith's pure E/B estimator in the flat-sky approximation.

    Parameters
    ----------
    nx : int
        Number of grid points in the x direction.
    ny : int
        Number of grid points in the y direction.
    D : array_like of float, shape (2,)
        Side lengths of the map in the x and y directions.
    QU : ndarray of float, shape (nx, ny, 2)
        Unmasked Q and U maps, with bounds ``(0:nx-1, 0:ny-1, 0:1)``.
    W : ndarray of float, shape (nx, ny)
        Window function, with bounds ``(0:nx-1, 0:ny-1)``.
    Wd : ndarray of float, shape (5, nx, ny), optional
        Precomputed derivatives of the window function:
        ``dW/dx``, ``dW/dy``, ``d^2W/dx^2``, ``d^2W/dxdy``, and
        ``d^2W/dy^2``. If not given, a zero array is used.

    Returns
    -------
    EB : ndarray of complex, shape (2, nx, ny)
        E and B modes on the two-dimensional Fourier grid, with bounds
        ``(0:1, 0:nx-1, 0:ny-1)``.

    Examples
    --------
    >>> EB = flatsky.ffttools.eb_separate(nx, ny, D, QU, W, Wd)
    """
    if Wd is None:
        Wd = np.zeros((5, nx, ny))
    return libflatsky.ffttools.eb_separate(nx, ny, D, QU, W, Wd)
    