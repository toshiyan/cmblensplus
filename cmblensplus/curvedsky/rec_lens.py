from ._core import lib_rec_lens


def qtt(lmax, rlmin, rlmax, fC, Tlm1, Tlm2, nside_t=0, gtype='', verbose=False, nthreads=0):
    """
    Reconstruct the CMB lensing potential and curl mode from the temperature quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing-potential alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        TT spectrum, with bounds ``0:rlmax``.
    Tlm1 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        First inverse-variance filtered temperature alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Tlm2 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Second inverse-variance filtered temperature alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        CMB lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode, or pseudo lensing-potential, alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_lens.qtt(
        lmax, rlmin, rlmax, fC, Tlm1, Tlm2,
        nside_t, gtype, verbose, nthreads=nthreads
    )


def qte(lmax, rlmin, rlmax, fC, Tlm, Elm, nside_t=0, gtype='', verbose=False, nthreads=0):
    """
    Reconstruct the CMB lensing potential and curl mode from the TE quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing-potential alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        TE spectrum, with bounds ``0:rlmax``.
    Tlm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered temperature alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Elm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered E-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        CMB lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode, or pseudo lensing-potential, alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_lens.qte(
        lmax, rlmin, rlmax, fC, Tlm, Elm,
        nside_t, gtype, verbose, nthreads=nthreads
    )


def qtb(lmax, rlmin, rlmax, fC, Tlm, Blm, nside_t=0, gtype='', verbose=False, nthreads=0):
    """
    Reconstruct the CMB lensing potential and curl mode from the TB quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing-potential alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        TE spectrum, with bounds ``0:rlmax``.
    Tlm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered temperature alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Blm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered B-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        CMB lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode, or pseudo lensing-potential, alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_lens.qtb(
        lmax, rlmin, rlmax, fC, Tlm, Blm,
        nside_t, gtype, verbose, nthreads=nthreads
    )


def qee(lmax, rlmin, rlmax, fC, Elm1, Elm2, nside_t=0, gtype='', verbose=False, nthreads=0):
    """
    Reconstruct the CMB lensing potential and curl mode from the EE quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing-potential alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        EE spectrum, with bounds ``0:rlmax``.
    Elm1 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        First inverse-variance filtered E-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Elm2 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Second inverse-variance filtered E-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        CMB lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode, or pseudo lensing-potential, alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_lens.qee(
        lmax, rlmin, rlmax, fC, Elm1, Elm2,
        nside_t, gtype, verbose, nthreads=nthreads
    )


def qeb(lmax, rlmin, rlmax, fC, Elm, Blm, nside_t=0, gtype='', verbose=False, nthreads=0):
    """
    Reconstruct the CMB lensing potential and curl mode from the EB quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing-potential alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        EE spectrum, with bounds ``0:rlmax``.
    Elm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered E-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Blm : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Inverse-variance filtered B-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        CMB lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode, or pseudo lensing-potential, alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_lens.qeb(
        lmax, rlmin, rlmax, fC, Elm, Blm,
        nside_t, gtype, verbose, nthreads=nthreads
    )


def qbb(lmax, rlmin, rlmax, fC, Blm1, Blm2, nside_t=0, gtype='', verbose=False, nthreads=0):
    """
    Reconstruct the CMB lensing potential and curl mode from the BB quadratic estimator.

    Parameters
    ----------
    lmax : int
        Maximum multipole of the output lensing-potential alms.
    rlmin : int
        Minimum CMB multipole used for reconstruction.
    rlmax : int
        Maximum CMB multipole used for reconstruction.
    fC : array_like of float, shape (rlmax + 1,)
        BB spectrum, with bounds ``0:rlmax``.
    Blm1 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        First inverse-variance filtered B-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    Blm2 : ndarray of complex, shape (rlmax + 1, rlmax + 1)
        Second inverse-variance filtered B-mode alm, with bounds
        ``(0:rlmax, 0:rlmax)``.
    nside_t : int, optional
        Nside for the convolution calculation. Default is 0.
    gtype : str, optional
        Type of output. Use ``'k'`` for convergence or ``''`` for lensing
        potential. Default is ``''``.
    verbose : bool, optional
        Whether to output messages. Default is False.
    nthreads : int, optional
        Number of threads. Default is 0.

    Returns
    -------
    glm : ndarray of complex, shape (lmax + 1, lmax + 1)
        CMB lensing-potential alm, with bounds ``(0:lmax, 0:lmax)``.
    clm : ndarray of complex, shape (lmax + 1, lmax + 1)
        Curl-mode, or pseudo lensing-potential, alm, with bounds
        ``(0:lmax, 0:lmax)``.
    """
    return lib_rec_lens.qbb(
        lmax, rlmin, rlmax, fC, Blm1, Blm2,
        nside_t, gtype, verbose, nthreads=nthreads
    )
    