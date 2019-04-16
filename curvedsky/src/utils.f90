!////////////////////////////////////////////////////!
! Utils for curvedsky analysis
!////////////////////////////////////////////////////!

module utils
  !from Healpix
  use alm_tools, only: alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  use fitstools, only: input_map, output_map
  use pix_tools, only: pix2ang_ring, ang2pix_ring, convert_ring2nest, convert_nest2ring
  use mask_tools, only: fill_holes_nest, dist2holes_nest
  !from F90/src_utils
  use random,    only: InitRandom, Gaussian1
  use constants, only: pi, iu
  use general,   only: str
  use pstool,    only: binned_ells, power_binning
  implicit none

  private alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  private input_map, output_map
  private pix2ang_ring, ang2pix_ring
  private initrandom, gaussian1
  private pi, iu
  private str
  private binned_ells, power_binning
  private fill_holes_nest, dist2holes_nest
  private convert_ring2nest, convert_nest2ring

contains


subroutine gauss1alm(lmax,cl,alm)
!*  Generating alm as a random Gaussian field whose power spectrum is cl.
!*  The output alm is given by a 2D array.
!*
!*  Args:
!*    - lmax (int)        : maximum multipole of the output alm
!*    - Cl[l] (double)    : angular power spectrum, with bounds (1:lmax)
!*
!*  Returns:
!*    - alm[l,m] (dcmplx) : random Gaussian alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(0:lmax) :: Cl
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l,m
  
  call initrandom(-1)
  alm = 0
  do l = 1, lmax
    alm(l,0) = Gaussian1()* sqrt(Cl(l))
    do m = 1, l
      alm(l,m) = cmplx(Gaussian1(),Gaussian1())*sqrt(Cl(l)/2)
    end do 
  end do

end subroutine gauss1alm


subroutine gauss2alm(lmax,cl1,cl2,xl,alm)
!*  Generating two alms as random Gaussian fields whose power spectra are cl1, cl2 and the cross spectrum is xl.
!*
!*  Args:
!*    - lmax (int)      : maximum multipole of the output alm
!*    - cl1[l] (double) : angular power spectrum of the 1st alm, with bounds (1:lmax)
!*    - cl2[l] (double) : angular power spectrum of the 2nd alm, with bounds (1:lmax)
!*    - xl[l] (double)  : cross-power spectrum between alm1 and alm2, with bounds (1:lmax)
!*
!*  Returns:
!*    - alm[2,l,m] (dcmplx): random Gaussian alms, with bounds (2,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(0:lmax) :: cl1, cl2, xl
  double complex, intent(out), dimension(2,0:lmax,0:lmax) :: alm
  !internal
  integer :: l, m
  double precision :: tamp, corr, xamp

  alm = 0
  call gauss1alm(lmax,cl1,alm(1,:,:))

  call initrandom(-1)
  do l = 2, lmax
    if (cl1(l)<=0d0)  stop 'error (gauss2alm): negative input cl1'
    corr = xl(l)/cl1(l)
    xamp = sqrt(cl2(l) - corr*xl(l))
    !m=0
    alm(2,l,0) = corr*alm(1,l,0) + Gaussian1()*xamp
    !m/=0
    do m =1, l
      alm(2,l,m) = corr*alm(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp/sqrt(2d0)
    end do
  end do

end subroutine gauss2alm


subroutine gaussTEB(lmax,TT,EE,BB,TE,alm)
!*  Generating T/E/B alms as random Gaussian fields whose power spectra are TT, EE, BB and the cross spectrum is TE.
!*
!*  Args:
!*    - lmax (int)     : maximum multipole of the output alms
!*    - TT[l] (double) : angular power spectrum of temperature, with bounds (1:lmax)
!*    - EE[l] (double) : angular power spectrum of E mode, with bounds (1:lmax)
!*    - BB[l] (double) : angular power spectrum of B mode, with bounds (1:lmax)
!*    - TE[l] (double) : TE cross-power spectrum, with bounds (1:lmax)
!*
!*  Returns:
!*    - alm[3,l,m] (dcmplx): random Gaussian T/E/B alms, with bounds (3,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(0:lmax) :: TT, EE, BB, TE
  double complex, intent(out), dimension(3,0:lmax,0:lmax) :: alm

  call gauss2alm(lmax,TT,EE,TE,alm(1:2,:,:))
  call gauss1alm(lmax,BB,alm(3,:,:))

end subroutine gaussTEB


subroutine gauss3alm(lmax,cl,alm)
!*  Generating three alms as random Gaussian fields whose covariance is given by cl[i,j].
!*
!*  Args:
!*    - lmax (int)         : maximum multipole of the output alm
!*    - cl[i,j,l] (double) : covariance between the gaussian fields, with bounds (3,3,1:lmax)
!*
!*  Returns:
!*    - alm[3,l,m] (dcmplx): random Gaussian alms, with bounds (3,0:lmax,0:lmax)
!*
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(3,3,lmax) :: cl
  double complex, intent(out), dimension(3,0:lmax,0:lmax) :: alm
  !integer
  integer :: l, m
  double precision :: sqrt2, A11, A21, A22, A31, A32, A33
  double complex :: g1, g2, g3

  sqrt2 = sqrt(2.)
  alm = 0d0
  call initrandom(-1)

  do l = 1, lmax
    !m=0
    A11 = 0d0; A21=0d0; A22=0d0; A31=0d0; A32=0d0; A33=0d0
    g1  = Gaussian1()
    g2  = Gaussian1()
    g3  = Gaussian1()
    A11 = dsqrt(cl(1,1,l))
    if (A11/=0) A21 = cl(1,2,l)/A11
    A22 = dsqrt(cl(2,2,l)-A21**2)
    if (A11/=0) A31 = cl(1,3,l)/A11
    if (A22/=0) A32 = (cl(2,3,l)-A31*A21)/A22
    A33 = dsqrt(cl(3,3,l)-A31**2-A32**2)
    alm(1,l,0) = g1*A11
    alm(2,l,0) = g1*A21 + g2*A22
    alm(3,l,0) = g1*A31 + g2*A32 + g3*A33
    A11 = 0d0;  A21=0d0; A22=0d0; A31=0d0; A32=0d0; A33=0d0
    !m/=0
    do m = 1, l
      g1  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g2  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g3  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      A11 = dsqrt(cl(1,1,l))
      if (A11/=0) A21 = cl(1,2,l)/A11
      A22 = dsqrt(cl(2,2,l)-A21**2)
      if (A11/=0) A31 = cl(1,3,l)/A11
      if (A22/=0) A32 = (cl(2,3,l)-A21*A31)/A22
      A33 = dsqrt(cl(2,2,l)-A31**2-A32**2)
      alm(1,l,m) = g1*A11
      alm(2,l,m) = g1*A21 + g2*A22
      alm(3,l,m) = g1*A31 + g2*A32 + g3*A33
    end do
  end do

end subroutine gauss3alm


subroutine gauss4alm(lmax,cl,alm)
!*  Generating four alms as random Gaussian fields whose covariance is given by cl[i,j].
!*
!*  Args:
!*    - lmax (int)         : maximum multipole of the output alm
!*    - cl[i,j,l] (double) : covariance between the gaussian fields, with bounds (4,4,1:lmax)
!*
!*  Returns:
!*    - alm[4,l,m] (dcmplx): random Gaussian alms, with bounds (4,0:lmax,0:lmax)
!*
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(4,4,lmax) :: cl
  double complex, intent(out), dimension(4,0:lmax,0:lmax) :: alm
  !integer
  integer :: l, m
  double precision :: sqrt2, A11, A21, A22, A31, A32, A33, A41, A42, A43, A44
  double complex :: g1, g2, g3, g4

  sqrt2 = sqrt(2.)
  alm = 0d0
  call initrandom(-1)

  do l = 1, lmax
    A11 = 0d0; A21=0d0; A22=0d0; A31=0d0; A32=0d0; A33=0d0; A41=0d0; A42=0d0; A43=0d0; A44=0d0
    g1  = Gaussian1()
    g2  = Gaussian1()
    g3  = Gaussian1()
    g4  = Gaussian1()
    A11 = dsqrt(cl(1,1,l))
    if (A11/=0) A21 = cl(1,2,l)/A11
    A22 = dsqrt(cl(2,2,l)-A21**2)
    if (A11/=0) A31 = cl(1,3,l)/A11
    if (A22/=0) A32 = (cl(2,3,l)-A31*A21)/A22
    A33 = dsqrt(cl(3,3,l)-A31**2-A32**2)
    if (A11/=0) A41 = cl(1,4,l)/A11
    if (A22/=0) A42 = (cl(2,4,l)-A41*A21)/A22
    if (A33/=0) A43 = (cl(3,4,l)-A41*A31-A42*A32)/A33
    A44 = dsqrt(cl(4,4,l)-A41**2-A42**2-A43**2)
    alm(1,l,0) = g1*A11
    alm(2,l,0) = g1*A21 + g2*A22
    alm(3,l,0) = g1*A31 + g2*A32 + g3*A33
    alm(4,l,0) = g1*A41 + g2*A42 + g3*A43 + g4*A44
    A11 = 0d0;  A21=0d0; A22=0d0; A31=0d0; A32=0d0; A33=0d0; A41=0d0; A42=0d0; A43=0d0; A44=0d0
    do m = 1, l
      g1  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g2  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g3  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g4  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      A11 = dsqrt(cl(1,1,l))
      if (A11/=0) A21 = cl(1,2,l)/A11
      A22 = dsqrt(cl(2,2,l)-A21**2)
      if (A11/=0) A31 = cl(1,3,l)/A11
      if (A22/=0) A32 = (cl(2,3,l)-A21*A31)/A22
      A33 = dsqrt(cl(2,2,l)-A31**2-A32**2)
      if (A11/=0) A41 = cl(1,4,l)/A11
      if (A22/=0) A42 = (cl(2,4,l)-A21*A41)/A22
      if (A33/=0) A43 = (cl(3,4,l)-A41*A31-A42*A32)/A33
      A44 = dsqrt(cl(4,4,l)-A41**2-A42**2-A43**2)
      alm(1,l,m) = g1*A11
      alm(2,l,m) = g1*A21 + g2*A22
      alm(3,l,m) = g1*A31 + g2*A32 + g3*A33
      alm(4,l,m) = g1*A41 + g2*A42 + g3*A43 + g4*A44
    end do
  end do

end subroutine gauss4alm


subroutine get_baseline(npix,nside_subpatch,QU,blmap)
!*  Calculate baseline of each subpatch. The subpatches have the same size.
!*  Written by Ryo Nagata.
!*
!*  Args:
!*    - npix (int)           : pixel number of the full map
!*    - nside_subpatch (int) : Nside of sub patch
!*    - QU[pix,2] (double)   : Q/U maps, with bounds (0:npix-1,2)
!*
!*  Returns:
!*    - blmap[pix,2] (double): baseline maps, with bounds (0:npix-1,2)
!*
  implicit none
  !I/O
  integer, intent(in) :: npix, nside_subpatch
  double precision, intent(in), dimension(0:npix-1,2) :: QU
  double precision, intent(out), dimension(0:npix-1,2) :: blmap
  !internal
  integer :: i, j, nside, npix_subpatch, npix_inpatch
  integer, allocatable, dimension(:) :: patch_map
  double precision, allocatable, dimension(:,:) :: patch_pmean
  double precision :: theta, phi

  nside = int(sqrt(npix/12d0))
  npix_subpatch = 12*nside_subpatch**2
  npix_inpatch  = npix/npix_subpatch

  allocate(patch_map(0:npix-1))
  allocate(patch_pmean(0:npix_subpatch-1,2))

  patch_pmean = 0d0

  ! compute baseline
  do i = 0, npix-1
    !identify subpatch
    call pix2ang_ring(nside, i, theta, phi)
    call ang2pix_ring(nside_subpatch, theta, phi, j)
    patch_map(i) = j
    !patch mean
    patch_pmean(j,:) = patch_pmean(j,:) + QU(i,:)/dble(npix_inpatch)
  end do

  ! obtain baseline at full pixel
  do i = 0, npix-1
    blmap(i,:) = patch_pmean(patch_map(i),:)
  end do

  ! baseline subtraction
  !Pmap = Pmap - blmap

  deallocate(patch_map,patch_pmean)

end subroutine get_baseline


subroutine get_winmap(nside_large,nside_small,ipix_pix,apod,win_out)
!*  Return apodization window for subpatch.
!*  Written by Ryo Nagata.
!*
!*  Args:
!*    - nside_large (int) : Nside of sub patch
!*    - nside_small (int) : full Nside
!*    - ipix_pix (int)    : pixel index of full map
!*    - apod (double)     : apodization length
!*
!*  Returns:
!*    - wind_out (double) : aporization window at ipix_pix 
!*
  implicit none
  !I/O
  integer, intent(in) :: nside_large, nside_small,ipix_pix
  double precision, intent(in)  :: apod
  double precision, intent(out) :: win_out
  !internal
  integer :: ipix_sub
  double precision :: dtheta, theta_pix, phi_pix, theta_sub, phi_sub

  call pix2ang_ring(nside_small, ipix_pix, theta_pix, phi_pix)
  call ang2pix_ring(nside_large, theta_pix, phi_pix, ipix_sub)
  call pix2ang_ring(nside_large, ipix_sub, theta_sub, phi_sub)

  dtheta = dacos(dcos(theta_pix)*dcos(theta_sub) &
      + dsin(theta_pix)*dcos(phi_pix)*dsin(theta_sub)*dcos(phi_sub)  &
      + dsin(theta_pix)*dsin(phi_pix)*dsin(theta_sub)*dsin(phi_sub) )
  call get_apod_window( dtheta, apod, win_out)

end subroutine get_winmap


subroutine get_apod_window(s,a,w)
!*  A sine apodization window
!*
!*  Args:
!*    - s (double) : distance from the center of the window
!*    - a (double) : apodization length, nothing (a=1) to all (a=0)
!*
!*  Returns:
!*    - w (double) : aporization window
!*
  implicit none
  !I/O
  double precision, intent(in) :: s, a
  double precision, intent(out) :: w
  !internal
  double precision :: s_in,s_out,s0

  s_out = 1d0
  s_in  = s_out*a
  if ( s < s_in  ) then 
    w=1d0
  else if ( s > s_out ) then
    w=0d0
  else
    w=(s_out-s)/(s_out-s_in)-dsin(2d0*pi*(s_out-s)/(s_out-s_in))/2d0/pi
  end if

end subroutine get_apod_window


subroutine cosin_healpix(npix,lmax,cosin) 
!*  cos(theta) as a function of the Healpix pixel index
!*
!*  Args:
!*    - npix (int) : pixel number of the desired map
!*    - lmax (int) : maximum multipole
!*
!*  Returns:
!*    - cosin[pix] (double) : cos(theta), with bounds (0:npix-1)
!*
  !I/O
  implicit none
  integer, intent(in) :: npix, lmax
  double precision, intent(out), dimension(0:npix-1) :: cosin
  !internal
  integer :: nside
  double complex :: alm(1,0:lmax,0:lmax)

  nside = int(sqrt(npix/12d0))
  alm = 0d0
  alm(1,1,0) = 1d0
  cosin = 0d0
  call alm2map(nside,lmax,lmax,alm,cosin)
  cosin = cosin*dsqrt(pi/0.75d0)

end subroutine cosin_healpix


subroutine eb_separate(npix,lmax,W,QUin,QUou)
!*  E/B mode seperation based on the chi-field estimator.
!*
!*  Args:
!*    - npix (int)           : pixel number of the desired map
!*    - lmax (int)           : maximum multipole used for the harmonic transform internally
!*    - W[pix,2] (double)    : window function satisfying zero 1st and 2nd derivatives at the boundary, with bounds (0:npix-1,2)
!*    - QUin[pix,2] (double) : input QU map, with bounds (0:npix-1,2)
!*
!*  Returns:
!*    - QUou[pix,2] (double) : E/B separated QU map, with bounds (0:npix-1,2)
!*
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax
  double precision, intent(in), dimension(0:npix-1,2) :: W, QUin
  double precision, intent(out), dimension(0:npix-1,2) :: QUou
  !internal
  integer :: i, l, nside
  double precision :: al, n1, n2, pd
  double precision, allocatable, dimension(:,:) :: W1,W2,P1,P0
  double complex, allocatable, dimension(:,:,:) :: alm1, alm2, wlm, tlm

  nside = int(sqrt(npix/12d0))

  !derivatives of window functions in harmonic space
  allocate(wlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,0,W,wlm)
  wlm = -wlm
  allocate(W1(0:npix-1,2),W2(0:npix-1,2),tlm(2,0:lmax,0:lmax))
  tlm = 0d0
  do l = 2, lmax
    al = dble(l)
    pd = dsqrt((al+1d0)*al)
    tlm(1,l,:) = - wlm(1,l,:)*pd
    tlm(2,l,:) = 0d0
  end do
  call alm2map_spin(nside,lmax,lmax,1,tlm,W1)
  do l = 2, lmax
    al = dble(l)
    pd = dsqrt((al+2d0)*(al**2-1d0)*al)
    tlm(1,l,:) = - wlm(1,l,:)*pd
    tlm(2,l,:) = 0d0
  end do
  call alm2map_spin(nside,lmax,lmax,2,tlm,W2)
  deallocate(tlm,wlm)

  !compute P1=conjg(W1)*(Q+/-iU), P0=conjg(W2)*(Q+/-iU)
  allocate(P1(0:npix-1,2),P0(0:npix-1,2))
  P1(:,1) = W1(:,1)*QUin(:,1)+W1(:,2)*QUin(:,2)
  P1(:,2) = W1(:,1)*QUin(:,2)-W1(:,2)*QUin(:,1)
  P0(:,1) = W2(:,1)*QUin(:,1)+W2(:,2)*QUin(:,2)
  P0(:,2) = W2(:,1)*QUin(:,2)-W2(:,2)*QUin(:,1)
  deallocate(W1,W2)

  allocate(alm1(2,0:lmax,0:lmax),alm2(2,0:lmax,0:lmax))
  alm1 = 0d0
  call map2alm_spin(nside,lmax,lmax,1,P1,alm1)
  call map2alm_spin(nside,lmax,lmax,0,P0,alm2)
  deallocate(P1,P0)

  do l = 2, lmax
    al = dble(l)
    n1 = 2d0/dsqrt((al+2d0)*(al-1d0))
    n2 = 1d0/dsqrt((al+2d0)*(al**2-1d0)*al)
    alm1(:,l,:) = n1*alm1(:,l,:) + n2*alm2(:,l,:)
  end do
  deallocate(alm2)

  allocate(P1(0:npix-1,2))
  call alm2map_spin(nside,lmax,lmax,2,alm1,P1)
  deallocate(alm1)

  ! P2 = conjg(W)*(Q+/-iU)
  QUou(:,1) = QUin(:,1)*W(:,1)
  QUou(:,2) = QUin(:,2)*W(:,1)

  ! total
  QUou = QUou + P1

end subroutine eb_separate


! From alm to angular power spectrum
subroutine alm2cl(lmax,cl,alm1,alm2,norm)
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(0:lmax) :: cl
  !optional
  double complex, intent(in), dimension(0:lmax,0:lmax), optional :: alm2
  double precision, intent(in), optional :: norm
  !f2py double complex :: alm2 = 0
  !f2py double precision :: norm = 1
  !docstr :: alm2 = alm1
  !intenral
  integer :: l

  cl = 0d0
  if (present(alm2).and..not.sum(abs(alm2))/=0) then
    do l = 1, lmax
      cl(l) = ( dble(alm1(l,0)*alm2(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm2(l,1:l))))/(2.*l+1.)
    end do
  else
    do l = 1, lmax
      cl(l) = ( dble(alm1(l,0)*alm1(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm1(l,1:l))))/(2.*l+1.)
    end do
  end if

  if (present(norm))  cl = cl*norm

end subroutine alm2cl


subroutine alm2bcl(bn,lmax,cb,alm1,alm2,oL,norm)
  implicit none
  integer, intent(in) :: bn, lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(bn) :: cb
  !optional
  integer, intent(in), optional, dimension(2) :: oL
  double precision, intent(in), optional :: norm
  double complex, intent(in), dimension(0:lmax,0:lmax), optional :: alm2
  !f2py integer :: oL = 0
  !f2py double precision norm = 1
  !f2py double complex :: alm2 = 0
  !docstr :: alm2 = alm1
  !internal
  integer :: npix, eL(2)
  double precision :: n
  double precision, allocatable :: Cl(:), bp(:), bc0(:)

  eL = (/1,lmax/)
  if (present(oL)) eL=oL

  n = 1d0
  if (present(norm)) n = norm

  allocate(Cl(lmax))
  if (present(alm2).and..not.sum(abs(alm2))/=0) then
    call alm2cl(lmax,cl,alm1,alm2,norm=n)
  else
    call alm2cl(lmax,cl,alm1,norm=n)
  end if

  allocate(bp(bn+1))
  call binned_ells(eL,bp)
  call power_binning(bp,eL,Cl,Cb)
  deallocate(bp,Cl)

end subroutine alm2bcl


subroutine apodize(npix,rmask,ascale,order,holeminsize,amask)
! partially use Healpix's process_mask code
!   ascale [deg] --- apodization length from the closest masked pixel
!   order        --- 1 for RING, otherwize NESTED
!   holeminsize [arcmin] --- minimum hole size (i.e., holes within this size in filled)
  implicit none
  integer, intent(in) :: npix
  double precision, intent(in) :: ascale
  double precision, intent(in), dimension(0:npix-1) :: rmask
  double precision, intent(out), dimension(0:npix-1) :: amask
  !optional
  integer, intent(in), optional :: order
  double precision, intent(in), optional :: holeminsize
  !f2py integer :: order = 1
  !f2py double precision :: holeminsize = 0
  !internal
  integer :: n, nside, hsize
  integer, allocatable :: mask(:)
  double precision :: x, y
  double precision, allocatable :: map(:,:)

  nside = int(sqrt(npix/12d0))

  allocate(map(0:npix-1,1)); map=0d0
  map(:,1) = rmask

  if (present(order).and.order == 1) then
     write(*,*) 'converting RING -> NESTED'
     call convert_ring2nest(nside, map)
  end if

  allocate(mask(0:npix-1))
  mask = NINT(map(:,1))

  !remove small holes
  hsize = nint( holeminsize / (4*pi/npix * (60d0*180d0/pi)**2) )
  if (hsize>0d0)  call fill_holes_nest(nside, hsize, mask, mask)

  !compute distance from valid to the closest invalid pixel
  write(*,*) 'compute distance'
  call dist2holes_nest(nside, mask, map(:,1))
  deallocate(mask)

  if (present(order).and.order == 1) then
     write(*,*) 'converting NESTED - > RING'
     call convert_nest2ring(nside, map)
  end if

  !compute apodization
  write(*,*) 'compute apodized window'
  amask = 0d0
  do n = 1, npix
    x = (1d0-dcos(map(n,1)))/(1d0-dcos(ascale*pi/180d0))
    y = min(1d0, dsqrt(x) )
    amask(n-1) = (y - dsin(2*pi*y)/(2*pi)) * rmask(n)
  end do

  deallocate(map)

end subroutine apodize


subroutine hp_map2alm(nside,lmax,mmax,map,alm)
  implicit none
  !I/O
  integer, intent(in) :: nside, lmax, mmax
  double precision, intent(in), dimension(:) :: map
  double complex, intent(out), dimension(0:lmax,0:mmax) :: alm
  !internal
  integer :: npix
  double complex :: tlm(1,0:lmax,0:mmax)

  npix = 12*nside**2

  if (size(map)/=npix) stop 'error (utils.hp_map2alm): size of map array is not 12*nside**2'

  tlm = 0d0
  call map2alm(nside,lmax,lmax,map,tlm)
  alm = tlm(1,:,:)

end subroutine hp_map2alm


subroutine hp_map2alm_spin(nside,lmax,mmax,spin,map0,map1,alm)
  implicit none
  !I/O
  integer, intent(in) :: nside, lmax, mmax, spin
  double precision, intent(in), dimension(:) :: map0, map1
  double complex, intent(out), dimension(2,0:lmax,0:mmax) :: alm
  !internal
  integer :: npix
  double precision :: S(0:12*nside**2-1,2)

  npix = 12*nside**2
  if (size(map0)/=npix) stop 'error (utils.hp_map2alm_spin): size of map0 array is not 12*nside**2'
  if (size(map1)/=npix) stop 'error (utils.hp_map2alm_spin): size of map1 array is not 12*nside**2'

  S(:,1) = map0
  S(:,2) = map1
  alm = 0d0
  call map2alm_spin(nside,lmax,lmax,spin,S,alm)

end subroutine hp_map2alm_spin


end module utils


