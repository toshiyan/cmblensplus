from ._core import lib_cninv


def cnfilter_freq(
    cl, bl, iNcov, maps, chn=1, lmaxs=[0], nsides=[0], itns=[1],
    eps=[1e-6], filter='W', inl=None, verbose=False, ro=50, stat=''
):
    """
    Optimally combine multiple-frequency CMB maps.

    The input maps must be beam-convolved. This routine deconvolves the beam
    during filtering, and returns filtered alms after beam deconvolution.

    The number of fields is denoted by ``n``: temperature only has ``n = 1``,
    polarization only has ``n = 2``, and temperature plus polarization has
    ``n = 3``. The number of frequency maps is denoted by ``mn``.

    Parameters
    ----------
    cl : array_like of float, shape (n, lmax + 1)
        Theory signal power spectra, with bounds ``(0:n-1, 0:lmax)``.
    bl : array_like of float, shape (mn, lmax + 1)
        Beam spectra, with bounds ``(0:mn-1, 0:lmax)``.
    iNcov : ndarray of float, shape (n, mn, npix)
        Inverse noise variance at each pixel, with bounds
        ``(0:n-1, 0:mn-1, 0:npix-1)``.
    maps : ndarray of float, shape (n, mn, npix)
        Beam-convolved T, Q, and U maps, with bounds
        ``(0:n-1, 0:mn-1, 0:npix-1)``.
    chn : int, optional
        Number of grids for the preconditioner. Use ``chn = 1`` for the
        diagonal preconditioner. Default is 1.
    lmaxs : array_like of int, optional
        Maximum multipoles at each preconditioning chain. ``lmaxs[0]`` is the
        input maximum multipole of ``cl``. Default is ``[0]``.
    nsides : array_like of int, optional
        Nsides of the preconditioner. ``nsides[0]`` should be consistent with
        the input map Nside. Default is ``[0]``.
    itns : array_like of int, optional
        Numbers of iterations for the preconditioning chains. Default is
        ``[1]``.
    eps : array_like of float, optional
        Convergence thresholds for the iterations. The iteration terminates
        when the residual fraction becomes smaller than ``eps``. Default is
        ``[1e-6]``.
    filter : {'', 'W'}, optional
        Filtering type. Use ``''`` for C-inverse filtering or ``'W'`` for
        Wiener filtering. Default is ``'W'``.
    inl : array_like of float, shape (n, mn, lmax + 1), optional
        Inverse noise spectra. If not given, the white-noise case is assumed.
    verbose : bool, optional
        Whether to output messages. Default is False.
    ro : int, optional
        Interval for printing the residual fraction. For example, ``ro = 2``
        outputs once every two iterations. Default is 50.
    stat : str, optional
        Realtime status filename containing the residual fraction. If empty,
        no status file is written. Default is ``''``.

    Returns
    -------
    xlm : ndarray of complex, shape (n, lmax + 1, lmax + 1)
        C-inverse or Wiener-filtered multipoles, with bounds
        ``(0:n-1, 0:lmax, 0:lmax)``.
    """
    if inl is None:
        lmax = len(cl[0]) - 1
        inl = 0 * iNcov[:, :, :lmax + 1]

    return lib_cninv.cnfilter_freq(
        cl, bl, iNcov, maps,
        chn=chn, lmaxs=lmaxs, nsides=nsides, itns=itns, eps=eps,
        filter=filter, inl=inl, verbose=verbose, ro=ro, stat=stat
    )


def cnfilter_kappa(
    cov, iNcov, maps, chn=1, lmaxs=[0], nsides=[0], itns=[1],
    eps=[1e-6], inl=None, verbose=False, ro=50, stat=''
):
    r"""
    Compute inverse-variance weighted multipoles for multiple kappa tracers.

    This routine computes :math:`(C + N)^{-1} x` for multiple mass-tracer
    kappa maps.

    Parameters
    ----------
    cov : array_like of float, shape (n, n, lmax + 1)
        Signal covariance matrix for each multipole, with bounds
        ``(0:n-1, 0:n-1, 0:lmax)``.
    iNcov : ndarray of float, shape (n, npix)
        Inverse noise variance at each pixel, with bounds
        ``(0:n-1, 0:npix-1)``.
    maps : ndarray of float, shape (n, npix)
        Input kappa maps, with bounds ``(0:n-1, 0:npix-1)``.
    chn : int, optional
        Number of grids for the preconditioner. Use ``chn = 1`` for the
        diagonal preconditioner. Default is 1.
    lmaxs : array_like of int, optional
        Maximum multipoles at each preconditioning chain. ``lmaxs[0]`` is the
        input maximum multipole of ``cov``. Default is ``[0]``.
    nsides : array_like of int, optional
        Nsides of the preconditioner. ``nsides[0]`` should be consistent with
        the input map Nside. Default is ``[0]``.
    itns : array_like of int, optional
        Numbers of iterations for the preconditioning chains. Default is
        ``[1]``.
    eps : array_like of float, optional
        Convergence thresholds for the iterations. The iteration terminates
        when the residual fraction becomes smaller than ``eps``. Default is
        ``[1e-6]``.
    inl : array_like of float, shape (n, lmax + 1), optional
        Inverse noise spectra for each mass map. If not given, the white-noise
        case is assumed.
    verbose : bool, optional
        Whether to output messages. Default is False.
    ro : int, optional
        Interval for printing the residual fraction. For example, ``ro = 2``
        outputs once every two iterations. Default is 50.
    stat : str, optional
        Realtime status filename containing the residual fraction. If empty,
        no status file is written. Default is ``''``.

    Returns
    -------
    xlm : ndarray of complex, shape (n, lmax + 1, lmax + 1)
        Wiener-filtered multipoles, with bounds ``(0:n-1, 0:lmax, 0:lmax)``.
    """
    if inl is None:
        lmax = len(cov[0, 0]) - 1
        inl = 0 * iNcov[:, :lmax + 1]

    return lib_cninv.cnfilter_kappa(
        cov, iNcov, maps,
        chn=chn, lmaxs=lmaxs, nsides=nsides, itns=itns, eps=eps,
        inl=inl, verbose=verbose, ro=ro, stat=stat
    )


def cnfilter_freq_nside(
    cl, bl0, bl1, iNcov0, iNcov1, maps0, maps1,
    chn=1, lmaxs=[0], nsides0=[0], nsides1=[0], itns=[1],
    eps=[1e-6], filter='W', inl=None, verbose=False, reducmn=0,
    ro=50, stat=''
):
    """
    Optimally combine multiple-frequency maps with two different Nsides.

    This routine is similar to :func:`cnfilter_freq`, but supports input maps
    with two different Nsides. The input maps must be beam-convolved. This
    routine deconvolves the beam during filtering, and returns filtered alms
    after beam deconvolution.

    Parameters
    ----------
    cl : array_like of float, shape (n, lmax + 1)
        Theory signal power spectra, with bounds ``(0:n-1, 0:lmax)``.
    bl0 : array_like of float, shape (mn, lmax + 1)
        Beam spectra for the first Nside set, with bounds
        ``(0:mn-1, 0:lmax)``.
    bl1 : array_like of float, shape (mn, lmax + 1)
        Beam spectra for the second Nside set, with bounds
        ``(0:mn-1, 0:lmax)``.
    iNcov0 : ndarray of float, shape (n, mn, npix0)
        Inverse noise variance for the first Nside set, with bounds
        ``(0:n-1, 0:mn-1, 0:npix0-1)``.
    iNcov1 : ndarray of float, shape (n, mn, npix1)
        Inverse noise variance for the second Nside set, with bounds
        ``(0:n-1, 0:mn-1, 0:npix1-1)``.
    maps0 : ndarray of float, shape (n, mn, npix0)
        Beam-convolved T, Q, and U maps for the first Nside set, with bounds
        ``(0:n-1, 0:mn-1, 0:npix0-1)``.
    maps1 : ndarray of float, shape (n, mn, npix1)
        Beam-convolved T, Q, and U maps for the second Nside set, with bounds
        ``(0:n-1, 0:mn-1, 0:npix1-1)``.
    chn : int, optional
        Number of grids for the preconditioner. Use ``chn = 1`` for the
        diagonal preconditioner. Default is 1.
    lmaxs : array_like of int, optional
        Maximum multipoles at each preconditioning chain. ``lmaxs[0]`` is the
        input maximum multipole of ``cl``. Default is ``[0]``.
    nsides0 : array_like of int, optional
        Nsides of the preconditioner for the first map set. ``nsides0[0]``
        should be consistent with the first input map Nside. Default is
        ``[0]``.
    nsides1 : array_like of int, optional
        Nsides of the preconditioner for the second map set. ``nsides1[0]``
        should be consistent with the second input map Nside. Default is
        ``[0]``.
    itns : array_like of int, optional
        Numbers of iterations for the preconditioning chains. Default is
        ``[1]``.
    eps : array_like of float, optional
        Convergence thresholds for the iterations. The iteration terminates
        when the residual fraction becomes smaller than ``eps``. Default is
        ``[1e-6]``.
    filter : {'', 'W'}, optional
        Filtering type. Use ``''`` for C-inverse filtering or ``'W'`` for
        Wiener filtering. Default is ``'W'``.
    inl : array_like of float, shape (n, mn, lmax + 1), optional
        Inverse noise spectra. If not given, the white-noise case is assumed.
    verbose : bool, optional
        Whether to output messages. Default is False.
    reducmn : int, optional
        Whether to reduce the number of maps per chain. Use 0 for no reduction,
        1 to combine maps with the same Nside inside the multigrid chain, or 2
        to additionally combine the two Nside map sets into a single map inside
        the second chain for ``chain >= 3``. Default is 0.
    ro : int, optional
        Interval for printing the residual fraction. For example, ``ro = 2``
        outputs once every two iterations. Default is 50.
    stat : str, optional
        Realtime status filename containing the residual fraction. If empty,
        no status file is written. Default is ``''``.

    Returns
    -------
    xlm : ndarray of complex, shape (n, lmax + 1, lmax + 1)
        C-inverse or Wiener-filtered multipoles, with bounds
        ``(0:n-1, 0:lmax, 0:lmax)``.
    """
    if inl is None:
        lmax = len(cl[0]) - 1
        inl = 0 * iNcov0[:, :, :lmax + 1]

    return lib_cninv.cnfilter_freq_nside(
        cl, bl0, bl1, iNcov0, iNcov1, maps0, maps1,
        chn=chn,
        lmaxs=lmaxs,
        nsides0=nsides0,
        nsides1=nsides1,
        itns=itns,
        eps=eps,
        filter=filter,
        inl=inl,
        verbose=verbose,
        reducmn=reducmn,
        ro=ro,
        stat=stat,
    )