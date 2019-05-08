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
!*  Generating alm as a random Gaussian field whose power spectrum is cl. The output alm is given by a 2D array.
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of the output alm
!*    :Cl [l] (double)    : Angular power spectrum, with bounds (0:lmax)
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Random Gaussian alm, with bounds (0:lmax,0:lmax)
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
!*    :lmax (int)       : Maximum multipole of the output alm
!*    :cl1 [l] (double) : Angular power spectrum of the 1st alm, with bounds (0:lmax)
!*    :cl2 [l] (double) : Angular power spectrum of the 2nd alm, with bounds (0:lmax)
!*    :xl [l] (double)  : Cross-power spectrum between alm1 and alm2, with bounds (0:lmax)
!*
!*  Returns:
!*    :alm [2,l,m] (dcmplx): Random Gaussian alms, with bounds (2,0:lmax,0:lmax)
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
!*    :lmax (int)      : Maximum multipole of the output alms
!*    :TT [l] (double) : Angular power spectrum of temperature, with bounds (0:lmax)
!*    :EE [l] (double) : Angular power spectrum of E mode, with bounds (0:lmax)
!*    :BB [l] (double) : Angular power spectrum of B mode, with bounds (0:lmax)
!*    :TE [l] (double) : TE cross-power spectrum, with bounds (0:lmax)
!*
!*  Returns:
!*    :alm [3,l,m] (dcmplx): Random Gaussian T/E/B alms, with bounds (3,0:lmax,0:lmax)
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
!*    :lmax (int)          : Maximum multipole of the output alm
!*    :cl [i,j,l] (double) : Covariance between the gaussian fields, with bounds (3,3,0:lmax)
!*
!*  Returns:
!*    :alm [3,l,m] (dcmplx): Random Gaussian alms, with bounds (3,0:lmax,0:lmax)
!*
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(3,3,0:lmax) :: cl
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
!*    :lmax (int)          : Maximum multipole of the output alm
!*    :cl [i,j,l] (double) : Covariance between the gaussian fields, with bounds (4,4,0:lmax)
!*
!*  Returns:
!*    :alm [4,l,m] (dcmplx): Random Gaussian alms, with bounds (4,0:lmax,0:lmax)
!*
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(4,4,0:lmax) :: cl
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
!*    :npix (int)            : pixel number of the full map
!*    :nside_subpatch (int)  : Nside of sub patch
!*    :QU [pix,2] (double)   : Q/U maps, with bounds (0:npix-1,2)
!*
!*  Returns:
!*    :blmap [pix,2] (double): baseline maps, with bounds (0:npix-1,2)
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
!*    :nside_large (int) : Nside of sub patch
!*    :nside_small (int) : full Nside
!*    :ipix_pix (int)    : pixel index of full map
!*    :apod (double)     : apodization length
!*
!*  Returns:
!*    :wind_out (double) : aporization window at ipix_pix 
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
!*    :s (double) : Distance from the center of the window
!*    :a (double) : Apodization length, nothing (a=1) to all (a=0)
!*
!*  Returns:
!*    :w (double) : Aporization window
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


subroutine eb_separate(npix,lmax,W,Q,U,Elm,Blm)
!*  E/B mode seperation based on the chi-field estimator. See e.g. Sec.III.2 of arXiv:1305.7441 for numerical implimentation.
!*
!*  Args:
!*    :npix (int)        : Pixel number of the desired map
!*    :lmax (int)        : Maximum multipole used for the harmonic transform internally
!*    :W[pix] (double)   : Window function satisfying zero 1st and 2nd derivatives at the boundary, with bounds (0:npix-1)
!*    :Q/U[pix] (double) : Input Q/U map, with bounds (0:npix-1)
!*
!*  Returns:
!*    :Elm/Blm[l,m]  (dcmplx) : Seperated E/B modes, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax
  double precision, intent(in), dimension(0:npix-1) :: W, Q, U
  double complex, intent(out), dimension(0:lmax,0:lmax) :: Elm, Blm
  !internal
  integer :: i, l, nside
  double precision :: al, n1(lmax), n2(lmax)
  double precision, allocatable, dimension(:,:) :: W1, W2, P2, P1, P0
  double complex, allocatable, dimension(:,:,:) :: alm0, alm1, alm2, wlm, tlm

  nside = int(sqrt(npix/12d0))

  !N1 and N2 defined as Ns = sqrt((l+s)!/(l-s)!)
  do l = 1, lmax
    al = dble(l)
    n1(l) = dsqrt((al+1d0)*al)
    n2(l) = dsqrt((al+2d0)*(al**2-1d0)*al)
  end do

  !//// derivatives of window functions in harmonic space ////!
  allocate(wlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,W,wlm)

  !compute del^1 W
  allocate(W1(0:npix-1,2),W2(0:npix-1,2),tlm(2,0:lmax,0:lmax))
  tlm = 0d0
  do l = 2, lmax
    tlm(1,l,:) = wlm(1,l,:)*n1(l)
  end do
  call alm2map_spin(nside,lmax,lmax,1,tlm,W1)

  !compute del^2 W
  tlm = 0d0
  do l = 2, lmax
    tlm(1,l,:) = wlm(1,l,:)*n2(l)
  end do
  call alm2map_spin(nside,lmax,lmax,2,tlm,W2)

  deallocate(tlm,wlm)

  !//// compute P1=conjg(W1)*(Q+iU), P0=conjg(W2)*(Q+iU) ////!
  allocate(P2(0:npix-1,2),P1(0:npix-1,2),P0(0:npix-1,2))
  P2(:,1) = W*Q
  P2(:,2) = W*U
  P1(:,1) = W1(:,1)*Q + W1(:,2)*U
  P1(:,2) = W1(:,1)*U - W1(:,2)*Q
  P0(:,1) = W2(:,1)*Q + W2(:,2)*U
  P0(:,2) = W2(:,1)*U - W2(:,2)*Q
  deallocate(W1,W2)

  ! P1^+ +/- iP^- = \sum_lm Y_lm^s alm^E +/- i alm^B
  allocate(alm0(2,0:lmax,0:lmax),alm1(2,0:lmax,0:lmax),alm2(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,2,P2,alm2)
  call map2alm_spin(nside,lmax,lmax,1,P1,alm1)
  call map2alm_spin(nside,lmax,lmax,0,P0,alm0)
  deallocate(P2,P1,P0)

  !compute purified E/B modes
  Elm = 0d0
  Blm = 0d0
  do l = 2, lmax
    Elm(l,:) = alm2(1,l,:) + (2*n1(l)/n2(l))*alm1(1,l,:) + (1d0/n2(l))*alm0(1,l,:)
    Blm(l,:) = alm2(2,l,:) + (2*n1(l)/n2(l))*alm1(2,l,:) + (1d0/n2(l))*alm0(2,l,:)
  end do
  deallocate(alm0,alm1,alm2)


end subroutine eb_separate


subroutine alm2cl(lmax,cl,alm1,alm2)
!*  From alm to angular power spectrum
!*
!*  Args:
!*    :lmax (int)           : Maximum multipole of the input alm
!*    :alm1 [l,m] (dcmplx)  : 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :alm2 [l,m] (dcmplx)  : 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1
!*
!*  Returns:
!*    :cl [l] (double) : Auto or cross angular power spectrum, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(0:lmax) :: cl
  !optional
  double complex, intent(in), optional, dimension(0:lmax,0:lmax) :: alm2
  !f2py double complex :: alm2 = 0
  !docstr :: alm2 = alm1

  !intenral
  integer :: l

  cl = 0d0
  if (present(alm2).and.sum(abs(alm2))/=0) then
    do l = 1, lmax
      cl(l) = ( dble(alm1(l,0)*alm2(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm2(l,1:l))))/(2.*l+1.)
    end do
  else
    do l = 1, lmax
      cl(l) = ( dble(alm1(l,0)*alm1(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm1(l,1:l))))/(2.*l+1.)
    end do
  end if

end subroutine alm2cl


subroutine alm2bcl(bn,lmax,cb,alm1,alm2,spc)
!*  From alm to angular power spectrum with multipole binning
!*
!*  Args:
!*    :bn (int)             : Number of multipole bins
!*    :lmax (int)           : maximum multipole of the input alm
!*    :alm1 [l,m] (dcmplx)  : 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :alm2 [l,m] (dcmplx)  : 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1
!*    :spc (str) : Specify bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing
!*
!*  Returns:
!*    :cb [bin] (double) : Auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)
!*
  implicit none
  integer, intent(in) :: bn, lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(bn) :: cb
  !optional
  character(4), intent(in), optional :: spc
  double complex, intent(in), optional, dimension(0:lmax,0:lmax) :: alm2
  !f2py character(4) :: spc = ''
  !f2py double complex :: alm2 = 0
  !docstr :: alm2 = alm1
  !internal
  integer :: npix, eL(2)
  double precision, allocatable :: Cl(:), bp(:), bc0(:)
  character(4) :: sp

  sp = ''
  if (present(spc)) sp = spc

  eL = (/1,lmax/)

  allocate(cl(0:lmax))
  if (present(alm2).and.sum(abs(alm2))/=0) then
    call alm2cl(lmax,cl,alm1,alm2)
  else
    call alm2cl(lmax,cl,alm1)
  end if

  allocate(bp(bn+1))
  call binned_ells(eL,bp,spc=sp)
  call power_binning(bp,eL,cl(1:lmax),Cb)
  deallocate(bp,Cl)

end subroutine alm2bcl


subroutine apodize(npix,rmask,ascale,order,holeminsize,amask)
!*  Compute apodized window function. Partially use Healpix's process_mask code.
!*
!*  Args:
!*   :npix (int)         : Number of pixel
!*   :rmask[pix] (double): Input window function, with bounds (0:pix-1)
!*   :ascale (double)    : Apodization length [deg] from the closest masked pixel
!*
!*  Args(optional):
!*   :order (int)         : Pixel order, 1 for RING (default), otherwize NESTED
!*   :holeminsize (double): Minimum hole size [arcmin] (i.e., holes within this size in filled), default to 0
!*
!*  Returns:
!*   :amask[pix] (double): Apodization window, with bounds (0:npix-1), using the same ordering as input
!*
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
!*  Ylm transform of the map to alm with the healpix (l,m) order
!*
!*  Args:
!*    :nside (int)         : Nside of the input map
!*    :lmax (int)          : Maximum multipole of the input alm
!*    :mmax (int)          : Maximum m of the input alm
!*    :map [pix] (double)  : Input map, with bounds (0:npix-1)
!*
!*  Returns:
!*    :alm1 [l,m] (dcmplx): Harmonic coefficient obtained from the input map, with bounds (0:lmax,0:mmax)
!*
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
!*  Spin Ylm transform of the map ( = map0 + i map1 ) to alm with the healpix (l,m) order. For example, if map0=Q, map1=U and spin=2, 
!*  the alm contains E-mode and B-mode. 
!*
!*  Args:
!*    :nside (int)         : Nside of the input map
!*    :lmax (int)          : Maximum multipole of the input alm
!*    :mmax (int)          : Maximum m of the input alm
!*    :spin (int)          : Spin of the transform
!*    :map0 [pix] (double) : Real part of the input map, with bounds (0:npix-1)
!*    :map1 [pix] (double) : Imaginary part of the input map, with bounds (0:npix-1)
!*
!*  Returns:
!*    :alm [2,l,m] (dcmplx): Parity-eve/odd harmonic coefficients obtained from the input map, with bounds (2,0:lmax,0:mmax)
!*
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


subroutine lm_healpy2healpix(lmpy,almpy,lmax,almpix)
!*  Transform healpy alm to healpix alm 
!*
!*  Args:
!*    :lmpy (int)           : Length of healpy alm
!*    :lmax (int)           : Maximum multipole
!*    :almpy[index] (dcmplx): Healpy alm, with bounds (0:lmpy-1)
!*
!*  Returns:
!*    :almpix [l,m] (dcmplx): Healpix alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  integer, intent(in) :: lmpy, lmax
  double complex, intent(in), dimension(0:lmpy-1) :: almpy
  double complex, intent(out), dimension(0:lmax,0:lmax) :: almpix
  integer :: l, m, i

  !i = m*(lmax+1) - m(m-1)/2 + (l-m)
  almpix = 0d0
  if (2*lmpy/=(lmax+1)*(lmax+2)) stop 'size of almpy is inconsistent with that of almpix'
  do m = 0, lmax
    do l = m, lmax
      i = m*lmax - int(m*(m-1d0)/2d0) + l
      almpix(l,m) = almpy(i)
    end do
  end do

end subroutine lm_healpy2healpix


subroutine cosin_healpix(npix,lmax,cosin)
!*  Return cos(theta) as a function of the Healpix pixel index
!*
!*  Args:
!*    :npix (int) : Pixel number of the desired map
!*    :lmax (int) : Maximum multipole
!*
!*  Returns:
!*    :cosin [pix] (double) : cosin(theta), with bounds (0:npix-1)
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


subroutine calc_mfs(bn,nu,lmax,walm,V,nside)
!*  Compute 2D Minkowski functionals
!*
!*  Args:
!*    :bn (int)            : Number of nu bins
!*    :nu [bin] (double)   : Nu bins, with bounds (bn)
!*    :lmax (int)          : Maximum multipole of the input walm
!*    :walm [l,m] (dcmplx) : Alm with filtering, possibly divided by the map variance, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional): 
!*    :nside (int): Nside of the intermediate map, default to lmax
!*
!*  Returns:
!*    :V [bin,type] (double): The three Minkowski functionals, V0, V1 and V2, at each nu bin, with bounds (bn,0:2)
!*
  implicit none 
  !I/O
  integer, intent(in) :: lmax, bn
  double precision, intent(in), dimension(bn) :: nu
  double complex, intent(in), dimension(0:lmax,0:lmax) :: walm
  double precision, intent(out), dimension(bn,0:2) :: V
  !optional
  integer, intent(in), optional :: nside
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !internal
  integer :: n, i, l, npix
  double precision :: dnu, dtt, dtp, dpp, atant, mcV(0:2)
  double precision, allocatable :: cost(:), der0(:), der1(:,:), der2(:,:)
  double complex, allocatable :: alm(:,:,:)

  npix = 12*nside**2
  dnu  = nu(2)-nu(1)

  ! compute derivatives 
  allocate(alm(1,0:lmax,0:lmax),der0(0:npix-1),der1(0:npix-1,2),der2(0:npix-1,3))
  alm(1,:,:) = walm
  call alm2map_der(nside,lmax,lmax,alm,der0,der1,der2)
  deallocate(alm)

  ! cosine
  allocate(cost(0:npix-1))
  call cosin_healpix(npix,lmax,cost) 

  !loop for nu
  do i = 1, bn
    mcV = 0d0
    do n = 0, npix-1
      !V0
      if (der0(n)>=nu(i))  mcV(0) = mcV(0) + 1d0
      if (abs(der0(n)-nu(i))<=dnu/2) then 
        !V1
        mcV(1) = mcV(1) + (0.25/dnu)*dsqrt(der1(n,1)**2+der1(n,2)**2)
        !V2
        atant  = cost(n)/dsqrt(1d0-cost(n)**2)
        dtt    = der2(n,1)
        dtp    = der2(n,2) - atant*der1(n,2)
        dpp    = der2(n,3) + atant*der1(n,1)
        mcV(2) = mcV(2) + (0.5/pi)*(1d0/dnu) * (2d0*der1(n,1)*der1(n,2)*dtt - der1(n,1)**2*dpp-der1(n,2)**2*dtt)/(der1(n,1)**2+der1(n,2)**2)
      end if
    end do
    V(i,0:2) = mcV(0:2)/dble(npix)
  end do

  deallocate(cost,der0,der1,der2)

end subroutine calc_mfs


end module utils


