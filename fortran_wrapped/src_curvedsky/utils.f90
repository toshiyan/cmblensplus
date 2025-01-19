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
  use random,    only: InitRandom, poisson, gaussian1
  use constants, only: pi, iu
  use general,   only: str, matrix_diag_tri
  use pstool,    only: binned_ells, power_binning
  use utilsgal,  only: nz_SF_scal, pz_SF_scal, gbias
  !from F90/src_hp
  !use hp_spht,   only: spht_map2alm
  implicit none

  private alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  private input_map, output_map
  private pix2ang_ring, ang2pix_ring
  private initrandom, gaussian1, poisson
  private pi, iu
  private str, matrix_diag_tri
  private binned_ells, power_binning
  private fill_holes_nest, dist2holes_nest
  private convert_ring2nest, convert_nest2ring
  private nz_SF_scal, pz_SF_scal, gbias
  !private spht_map2alm

contains


subroutine gauss1alm(lmax,cl,alm,lmin)
!*  Generating alm as a random Gaussian field whose power spectrum is cl. The output alm is given by a 2D array.
!*
!*  Args:
!*    :cl [l] (double) : Angular power spectrum, with bounds (0:lmax)
!*
!*  Args(optional):
!*    :lmin (int)      : Minimum multipole of output alm (default: 2)
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Random Gaussian alm, with bounds (0:lmax,0:lmax)
!*
  !f2py intent(in) lmax, lmin, cl
  !f2py intent(out) alm
  !f2py depend(lmax) cl, alm
  implicit none
  !I/O
  integer, intent(in) :: lmax, lmin
  double precision, intent(in), dimension(0:lmax) :: cl
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, m
  !rmargs :: lmax
  !opt4py :: lmin = 2
  !add2py :: lmax = len(cl) - 1
  !add2py :: if cl[lmin] <= 0: raise SystemExit('ERROR in "gauss1alm": cl[lmin] should be positive')
  
  call initrandom(-1)

  alm = 0
  do l = lmin, lmax
    alm(l,0) = gaussian1()* dsqrt(cl(l))
    do m = 1, l
      alm(l,m) = cmplx(gaussian1(),gaussian1())*dsqrt(cl(l)/2d0)
    end do 
  end do

end subroutine gauss1alm

subroutine gauss2alm(lmax,cl1,cl2,xl,alm,lmin,flm)
!*  Generating two alms as random Gaussian fields whose power spectra are cl1, cl2 and the cross spectrum is xl.
!*
!*  Args:
!*    :cl1 [l] (double) : Angular power spectrum of the 1st alm, with bounds (0:lmax)
!*    :cl2 [l] (double) : Angular power spectrum of the 2nd alm, with bounds (0:lmax)
!*    :xl [l] (double)  : Cross-power spectrum between alm1 and alm2, with bounds (0:lmax)
!*
!*  Args(optional):
!*    :lmin (int)         : Minimum multipole of output alm (default: 2)
!*    :flm [l,m] (dcmplx) : Constrained realiation of alms whose power spectrum is cl1
!*
!*  Returns:
!*    :alm [2,l,m] (dcmplx): Random Gaussian alms, with bounds (2,0:lmax,0:lmax)
!*
  !f2py intent(in) lmax, lmin, cl1, cl2, xl, flm
  !f2py intent(out) alm
  !f2py depend(lmax) cl1, cl2, xl, flm, alm
  implicit none
  !I/O
  integer, intent(in) :: lmax, lmin
  double precision, intent(in), dimension(0:lmax) :: cl1, cl2, xl
  double complex, intent(in), dimension(0:lmax,0:lmax) :: flm
  double complex, intent(out), dimension(2,0:lmax,0:lmax) :: alm
  !internal
  integer :: l, m
  double precision :: tamp, corr, xamp
  double complex, allocatable :: talm(:,:)
  !rmargs :: lmax
  !opt4py :: lmin = 2
  !opt4py :: flm = None
  !add2py :: if flm is None: flm = numpy.zeros((lmax+1,lmax+1),dtype=numpy.complex128)
  !add2py :: lmax = len(cl1) - 1
  !add2py :: if numpy.any(xl[lmin:]**2/cl1[lmin:]/cl2[lmin:]>1): raise SystemExit('ERROR in "gauss2alm": correlation coefficient is larger than 1')

  !initialize
  alm = 0

  if (sum(abs(flm)) == 0d0) then
    !generate alm from cl1
    allocate(talm(0:lmax,0:lmax))
    call gauss1alm(lmax,cl1,talm,lmin)
    alm(1,:,:) = talm
    deallocate(talm)
  else
    alm(1,:,:) = flm
  end if

  !generate correlated alm
  call initrandom(-1)
  do l = lmin, lmax
    if (cl1(l)<=0d0) then
        write(*,*) 'error (gauss2alm): negative input cl1'
        stop
    end if
    corr = xl(l)/cl1(l)
    if (cl2(l)<corr*xl(l)) then
        write(*,*) 'error (gauss2alm): correlation coefficient > 1 at specific multipole'
        stop
    end if
    xamp = sqrt(cl2(l) - corr*xl(l))
    !m=0
    alm(2,l,0) = corr*alm(1,l,0) + Gaussian1()*xamp
    !m/=0
    do m = 1, l
      alm(2,l,m) = corr*alm(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp/sqrt(2d0)
    end do
  end do

end subroutine gauss2alm

subroutine gaussTEB(lmax,TT,EE,BB,TE,alm,lmin)
!*  Generating T/E/B alms as random Gaussian fields whose power spectra are TT, EE, BB and the cross spectrum is TE.
!*
!*  Args:
!*    :TT [l] (double) : Angular power spectrum of temperature, with bounds (0:lmax)
!*    :EE [l] (double) : Angular power spectrum of E mode, with bounds (0:lmax)
!*    :BB [l] (double) : Angular power spectrum of B mode, with bounds (0:lmax)
!*    :TE [l] (double) : TE cross-power spectrum, with bounds (0:lmax)
!*
!*  Args(optional):
!*    :lmin (int)      : Minimum multipole of output alm (default: 2)
!*
!*  Returns:
!*    :alm [3,l,m] (dcmplx): Random Gaussian T/E/B alms, with bounds (3,0:lmax,0:lmax)
!*
  !f2py intent(in) lmax, lmin, TT, EE, BB, TE
  !f2py intent(out) alm
  !f2py depend(lmax) TT, EE, BB, TE, alm
  implicit none
  !I/O
  integer, intent(in) :: lmax, lmin
  double precision, intent(in), dimension(0:lmax) :: TT, EE, BB, TE
  double complex, intent(out), dimension(3,0:lmax,0:lmax) :: alm
  !internal
  double complex, allocatable, dimension(:,:) :: almb
  double complex, allocatable, dimension(:,:,:) :: almte
  !rmargs :: lmax
  !opt4py :: lmin = 2
  !add2py :: lmax = len(TT) - 1
  !add2py :: if numpy.any(TT[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": TT[lmin:] should be positive')
  !add2py :: if numpy.any(EE[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": EE[lmin:] should be positive')
  !add2py :: if numpy.any(BB[lmin:]<=0): raise SystemExit('ERROR in "gaussTEB": BB[lmin:] should be positive')
  !add2py :: if numpy.any(TE[2:]**2/TT[2:]/EE[2:]>1): raise SystemExit('ERROR in "gaussTEB": correlation coefficient is larger than 1')
  
  allocate(almte(2,0:lmax,0:lmax),almb(0:lmax,0:lmax))
  almb = 0d0
  call gauss2alm(lmax,TT,EE,TE,almte,lmin,almb)
  alm(1:2,:,:) = almte
  deallocate(almte)
  
  call gauss1alm(lmax,BB,almb,2) !lmin=2 for polarization
  alm(3,:,:) = almb
  deallocate(almb)

end subroutine gaussTEB

subroutine gaussalm(n,lmax,cl,alm,ilm)
!*  Generating alms as random Gaussian fields whose covariance is given by cl[i,j].
!*
!*  Args:
!*    :cl [i,j,l] (double) : Covariance between the gaussian fields, with bounds (n,n,0:lmax)
!*
!*  Args(optional):
!*    :ilm [l,m] (dcmplx)  : Input alm for the cl[0,0] element (default to None). The other alms are generated to be correlated with ilm. 
!*
!*  Returns:
!*    :alm [i,l,m] (dcmplx): Random Gaussian alms, with bounds (n,0:lmax,0:lmax)
!*
  !f2py intent(in) n, lmax, cl, ilm
  !f2py intent(out) alm
  !f2py depend(n) cl, alm
  !f2py depend(lmax) cl, ilm, alm
  implicit none
  !I/O
  integer, intent(in) :: n, lmax
  double precision, intent(in), dimension(n,n,0:lmax) :: cl
  double complex, intent(in), dimension(0:lmax,0:lmax) :: ilm
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: alm
  !integer
  integer :: l, m, i
  double precision :: prod, sumilm
  double precision, dimension(n,n) :: A
  double complex, dimension(n) :: g
  !rmargs :: n, lmax
  !opt4py :: n = None
  !opt4py :: lmax = None
  !opt4py :: ilm = None
  !add2py :: if n is None:    n    = len(cl[:,0,0])
  !add2py :: if lmax is None: lmax = len(cl[0,0,:]) - 1
  !add2py :: if ilm is None:  ilm  = [[0 for x in range(lmax+1)] for y in range(lmax+1)] 

  alm = 0d0
  call initrandom(-1)

  sumilm = sum(ilm)

  do l = 1, lmax

    ! check diagonal
    prod = 1d0
    do i = 1, n
      prod = prod * cl(i,i,l)
    end do
    if (prod==0d0) cycle

    call matrix_diag_tri(n,cl(:,:,l),A)

    do i = 1, n
      g(i) = Gaussian1()
      if (sumilm/=0 .and. i==1 .and. A(1,1)/=0d0)  g(1) = ilm(l,0)/A(1,1) ! replace to ilm
      alm(i,l,0) = sum(g(1:i)*A(i,1:i))
    end do

    do m = 1, l

      do i = 1, n
        g(i) = cmplx(Gaussian1(),Gaussian1())/dsqrt(2d0)
        if (sumilm/=0 .and. i==1 .and. A(1,1)/=0d0)  g(1) = ilm(l,m)/A(1,1) ! replace to ilm
        alm(i,l,m) = sum(g(1:i)*A(i,1:i))
      end do

    end do
  end do

end subroutine gaussalm

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
  !f2py intent(in) npix, nside_subpatch, QU
  !f2py intent(out) blmap
  !f2py depend(npix) QU, blmap
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
  !f2py intent(in) nside_large, nside_small, ipix_pix, apod
  !f2py intent(out) win_out
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
  !f2py intent(in) s, a
  !f2py intent(out) w
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
!*    :lmax (int)        : Maximum multipole used for the harmonic transform internally
!*    :W[pix] (double)   : Window function satisfying zero 1st and 2nd derivatives at the boundary, with bounds (0:npix-1)
!*    :Q/U[pix] (double) : Input Q/U map already multiplied by W, with bounds (0:npix-1)
!*
!*  Returns:
!*    :Elm/Blm[l,m]  (dcmplx) : Seperated E/B modes, with bounds (0:lmax,0:lmax)
!*
  !f2py intent(in) npix, lmax, W, Q, U
  !f2py intent(out) Elm, Blm
  !f2py depend(npix) W, Q, U
  !f2py depend(lmax) Elm, Blm
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax
  double precision, intent(in), dimension(0:npix-1) :: W, Q, U
  double complex, intent(out), dimension(0:lmax,0:lmax) :: Elm, Blm
  !internal
  integer :: i, n, l, nside
  double precision :: al, n1(lmax), n2(lmax)
  double precision, allocatable, dimension(:,:) :: W1, W2, P2, P1, P0
  double complex, allocatable, dimension(:,:,:) :: alm0, alm1, alm2, wlm, tlm
  !rmargs :: npix
  !add2py :: npix = len(W)

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

  allocate(W1(0:npix-1,2),W2(0:npix-1,2),tlm(2,0:lmax,0:lmax))
  !compute del^1 W
  tlm = 0d0
  do l = 1, lmax
    tlm(1,l,0:l) = wlm(1,l,0:l)*n1(l)
  end do
  call alm2map_spin(nside,lmax,lmax,1,tlm,W1)

  !compute del^2 W
  tlm = 0d0
  do l = 2, lmax
    tlm(1,l,0:l) = wlm(1,l,0:l)*n2(l)
  end do
  call alm2map_spin(nside,lmax,lmax,2,tlm,W2)
  deallocate(tlm,wlm)

  !inverse mask
  do n = 0, npix-1
    if (W(n)==0d0) then
      W1(n,:) = 0d0
      W2(n,:) = 0d0
    else
      W1(n,:) = W1(n,:)/W(n)
      W2(n,:) = W2(n,:)/W(n)
    end if
  end do

  !//// compute P1=conjg(W1)*(Q+iU), P0=conjg(W2)*(Q+iU) ////!
  allocate(P2(0:npix-1,2),P1(0:npix-1,2),P0(0:npix-1,2))
  P2(:,1) = Q
  P2(:,2) = U
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
!*    :alm1 [l,m] (dcmplx)  : 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :alm2 [l,m] (dcmplx)  : 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1
!*
!*  Returns:
!*    :cl [l] (double) : Auto or cross angular power spectrum, with bounds (0:lmax)
!*
  !f2py intent(in) lmax, alm1, alm2
  !f2py intent(out) cl
  !f2py depend(lmax) alm1, cl, alm2
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(0:lmax) :: cl
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm2
  !intenral
  integer :: l
  !rmargs :: lmax
  !opt4py :: alm2 = None
  !add2py :: lmax = len(alm1[:,0]) - 1
  !add2py :: if alm2 is None:  alm2 = alm1

  cl = 0d0
  do l = 1, lmax
    cl(l) = ( dble(alm1(l,0)*alm2(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm2(l,1:l))))/(2.*l+1.)
  end do

end subroutine alm2cl

subroutine alm2bcl(bn,lmax,cb,alm1,alm2,spc)
!*  From alm to angular power spectrum with multipole binning
!*
!*  Args:
!*    :bn (int)             : Number of multipole bins
!*    :alm1 [l,m] (dcmplx)  : 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :alm2 [l,m] (dcmplx)  : 2nd harmonic coefficient, with bounds (0:lmax,0:lmax), default to alm1
!*    :spc (str) : Specify bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing
!*
!*  Returns:
!*    :cb [bin] (double) : Auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)
!*
  !f2py intent(in) bn, lmax, alm1, spc, alm2
  !f2py intent(out) cb
  !f2py depend(lmax) alm1, alm2
  !f2py depend(bn) cb
  implicit none
  integer, intent(in) :: bn, lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(bn) :: cb
  character(4), intent(in) :: spc
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm2
  !internal
  integer :: npix, eL(2)
  double precision, allocatable :: Cl(:), bp(:)
  !rmargs :: lmax
  !opt4py :: spc = ''
  !opt4py :: alm2 = None
  !add2py :: lmax = len(alm1[:,0]) - 1
  !add2py :: if alm2 is None:  alm2 = alm1

  eL = (/1,lmax/)

  allocate(cl(0:lmax))
  call alm2cl(lmax,cl,alm1,alm2)

  allocate(bp(bn+1))
  call binned_ells(eL,bp,spc=spc)
  call power_binning(bp,eL,cl(1:lmax),Cb)
  deallocate(bp,Cl)

end subroutine alm2bcl

subroutine alm2rho(lmax,rho,alm1,alm2)
!*  Compute correlation coefficients between two alms
!*
!*  Args:
!*    :alm1 [l,m] (dcmplx)  : 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
!*    :alm2 [l,m] (dcmplx)  : 2nd harmonic coefficient, with bounds (0:lmax,0:lmax)
!*
!*  Returns:
!*    :rho [l] (double) : Auto or cross angular power spectrum, with bounds (0:lmax)
!*
  !f2py intent(in) lmax, alm1, alm2
  !f2py intent(out) rho
  !f2py depend(lmax) alm1, alm2, rho
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm2
  double precision, intent(out), dimension(0:lmax) :: rho
  !intenral
  integer :: l
  double precision :: cl(3,0:lmax)
  !rmargs :: lmax
  !add2py :: lmax = len(alm1[:,0]) - 1

  call alm2cl(lmax,cl(1,:),alm1,alm1)
  call alm2cl(lmax,cl(2,:),alm2,alm2)
  call alm2cl(lmax,cl(3,:),alm1,alm2)
  
  rho = 0d0
  do l = 0, lmax
    if (cl(1,l)*cl(2,l)==0d0)  cycle
    rho(l) = cl(3,l)**2/(cl(1,l)*cl(2,l))
  end do

end subroutine alm2rho

subroutine alm2cov(n,lmax,cov,alm)
!*  Compute correlation coefficients between two alms
!*
!*  Args:
!*    :alm [n,l,m] (dcmplx) : 1st harmonic coefficient, with bounds (0:lmax,0:lmax)
!*
!*  Returns:
!*    :cov [n,n,l] (double) : Auto and cross angular power spectra between alm[i] and alm[j]
!*
  !f2py intent(in) lmax, n, alm
  !f2py intent(out) cov
  !f2py depend(lmax) alm, cov
  !f2py depend(n) alm, cov
  implicit none
  !I/O
  integer, intent(in) :: lmax, n
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: alm
  double precision, intent(out), dimension(n,n,0:lmax) :: cov
  !intenral
  integer :: i1, i2
  !rmargs :: n, lmax
  !add2py :: n = len(alm[:,0,0])
  !add2py :: lmax = len(alm[0,:,0]) - 1

  do i1 = 1, n
    do i2 = i1, n
      call alm2cl(lmax,cov(i1,i2,:),alm(i1,:,:),alm(i2,:,:))
      cov(i2,i1,:) = cov(i1,i2,:)
    end do
  end do

end subroutine alm2cov

subroutine apodize(npix,rmask,ascale,order,holeminsize,amask)
!*  Compute apodized window function. Partially use Healpix's process_mask code.
!*
!*  Args:
!*    :rmask[pix] (double): Input window function, with bounds (0:pix-1). Pixels at rmask=0 is considered as masked pixels. 
!*    :ascale (double)    : Apodization length [deg] from the closest masked pixel
!*
!*  Args(optional):
!*    :order (int)         : Pixel order, 1 for RING (default), otherwize NESTED
!*    :holeminsize (double): Minimum hole size [arcmin] (i.e., holes within this size is filled), default to 0
!*
!*  Returns:
!*    :amask[pix] (double): Apodization window, with bounds (0:npix-1), using the same ordering as input
!*
  !f2py intent(in) npix, ascale, rmask, order, holeminsize
  !f2py intent(out) amask
  !f2py depend(npix) rmask, amask
  implicit none
  !I/O
  integer, intent(in) :: npix
  double precision, intent(in) :: ascale
  double precision, intent(in), dimension(0:npix-1) :: rmask
  double precision, intent(out), dimension(0:npix-1) :: amask
  integer, intent(in) :: order
  double precision, intent(in) :: holeminsize
  !internal
  integer :: n, nside, hsize
  integer, allocatable :: mask(:)
  double precision :: x, y
  double precision, allocatable :: map(:,:), rmap(:)
  !rmargs :: npix
  !opt4py :: order = 1
  !opt4py :: holeminsize = 0
  !add2py :: npix = len(rmask)

  nside = int(sqrt(npix/12d0))

  allocate(map(0:npix-1,1),rmap(0:npix-1)); map=0d0; rmap=0d0
  map(:,1) = rmask

  if (order == 1) then
     write(*,*) 'converting RING -> NESTED'
     call convert_ring2nest(nside, map)
  end if

  allocate(mask(0:npix-1)); mask=0
  do n = 0, npix-1
    if (map(n,1)/=0d0)  mask(n) = int(map(n,1)/map(n,1))
  end do

  !remove small holes
  hsize = nint( holeminsize / (4*pi/npix * (60d0*180d0/pi)**2) )
  if (hsize>0d0)  call fill_holes_nest(nside, hsize, mask, mask)

  !compute distance from valid to the closest invalid pixel
  write(*,*) 'compute distance'
  call dist2holes_nest(nside, mask, rmap)
  map(:,1) = rmap
  deallocate(mask)

  if (order == 1) then
     write(*,*) 'converting NESTED - > RING'
     call convert_nest2ring(nside, map)
  end if

  !compute apodization
  write(*,*) 'compute apodized window'
  amask = 0d0
  do n = 0, npix-1
    x = (1d0-dcos(map(n,1)))/(1d0-dcos(ascale*pi/180d0))
    y = min(1d0, dsqrt(x) )
    amask(n) = (y - dsin(2*pi*y)/(2*pi)) * rmask(n)
  end do

  deallocate(map)

end subroutine apodize

subroutine hp_alm2map(npix,lmax,mmax,alm,map)
!*  Ylm transform of the map to alm with the healpix (l,m) order
!*
!*  Args:
!*    :nside (int)        : Nside of the input map
!*    :alm [l,m] (dcmplx) : Harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
!*
!*  Returns:
!*    :map [pix] (double) : Transformed map, with bounds (0:npix-1)
!*
  !f2py intent(in) npix, lmax, mmax, alm
  !f2py intent(out) map
  !f2py depend(lmax) alm
  !f2py depend(mmax) alm
  !f2py depend(npix) map
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, mmax
  double complex, intent(in), dimension(0:lmax,0:mmax) :: alm
  double precision, intent(out), dimension(0:npix-1) :: map
  !internal
  integer :: nside
  double complex, allocatable :: tlm(:,:,:)
  double precision, allocatable :: tmap(:)
  !chargs :: npix -> nside
  !rmargs :: lmax, mmax
  !add2py :: lmax = len(alm[:,0]) - 1
  !add2py :: mmax = len(alm[0,:]) - 1
  !add2py :: npix = 12*nside**2

  nside = int(dsqrt(npix/12d0))

  allocate(tlm(1,0:lmax,0:lmax))
  tlm(1,:,:) = alm
  call alm2map(nside,lmax,mmax,tlm,map)
  deallocate(tlm)

end subroutine hp_alm2map

subroutine hp_alm2map_spin(npix,lmax,mmax,spin,elm,blm,map0,map1)
!*  Ylm transform of the map to alm with the healpix (l,m) order
!*
!*  Args:
!*    :nside (int)         : Nside of the input map
!*    :spin (int)         : Spin of the transform
!*    :elm [l,m] (dcmplx) : Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
!*    :blm [l,m] (dcmplx) : Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
!*
!*  Returns:
!*    :map0 [pix] (double): Real part of the transformed map (Q-like map), with bounds (0:npix-1)
!*    :map1 [pix] (double): Imaginary part of the transformed map (U-like map), with bounds (0:npix-1)
!*
  !f2py intent(in) npix, lmax, mmax, spin, elm, blm
  !f2py intent(out) map0, map1
  !f2py depend(lmax) elm, blm
  !f2py depend(mmax) elm, blm
  !f2py depend(npix) map0, map1
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, mmax, spin
  double complex, intent(in), dimension(0:lmax,0:mmax) :: elm, blm
  double precision, intent(out), dimension(0:npix-1) :: map0, map1
  !internal
  integer :: nside
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: alm(:,:,:)
  !chargs :: npix -> nside
  !rmargs :: lmax, mmax
  !add2py :: lmax = len(elm[:,0]) - 1
  !add2py :: mmax = len(elm[0,:]) - 1
  !add2py :: npix = 12*nside**2

  nside = int(dsqrt(npix/12d0))

  allocate(alm(2,0:lmax,0:mmax),map(0:npix-1,2))
  alm(1,:,:) = elm
  alm(2,:,:) = blm
  call alm2map_spin(nside,lmax,mmax,spin,alm,map)
  map0 = map(:,1)
  map1 = map(:,2)
  deallocate(alm,map)

end subroutine hp_alm2map_spin

subroutine hp_map2alm(npix,lmax,mmax,map,alm)
!*  Ylm transform of the map to alm with the healpix (l,m) order
!*
!*  Args:
!*    :lmax (int)          : Maximum multipole of the input alm
!*    :mmax (int)          : Maximum m of the input alm
!*    :map [pix] (double)  : Input map, with bounds (0:npix-1)
!*
!*  Returns:
!*    :alm [l,m] (dcmplx): Harmonic coefficient obtained from the input map, with bounds (0:lmax,0:mmax)
!*
  !f2py intent(in) npix, lmax, mmax, map
  !f2py intent(out) alm
  !f2py depend(npix) map
  !f2py depend(lmax) alm
  !f2py depend(mmax) alm
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, mmax
  double precision, intent(in), dimension(0:npix-1) :: map
  double complex, intent(out), dimension(0:lmax,0:mmax) :: alm
  !internal
  integer :: nside
  double complex, allocatable :: tlm(:,:,:)
  !rmargs :: npix
  !add2py :: npix = len(map)

  nside = int(dsqrt(npix/12d0))

  allocate(tlm(1,0:lmax,0:mmax)); tlm = 0d0
  call map2alm(nside,lmax,mmax,map,tlm)
  alm = tlm(1,:,:)
  deallocate(tlm)

end subroutine hp_map2alm

subroutine hp_map2alm_spin(npix,lmax,mmax,spin,map0,map1,alm)
!*  Spin Ylm transform of the map ( = map0 + i map1 ) to alm with the healpix (l,m) order. For example, if map0=Q, map1=U and spin=2, 
!*  the alm contains E-mode and B-mode. 
!*
!*  Args:
!*    :lmax (int)          : Maximum multipole of the input alm
!*    :mmax (int)          : Maximum m of the input alm
!*    :spin (int)          : Spin of the transform
!*    :map0 [pix] (double) : Real part of the input map, with bounds (0:npix-1)
!*    :map1 [pix] (double) : Imaginary part of the input map, with bounds (0:npix-1)
!*
!*  Returns:
!*    :alm [2,l,m] (dcmplx): Parity-eve/odd harmonic coefficients obtained from the input map, with bounds (2,0:lmax,0:mmax)
!*
  !f2py intent(in) npix, lmax, mmax, spin, map0, map1
  !f2py intent(out) alm
  !f2py depend(npix) map0, map1
  !f2py depend(lmax) alm
  !f2py depend(mmax) alm
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, mmax, spin
  double precision, intent(in), dimension(0:npix-1) :: map0, map1
  double complex, intent(out), dimension(2,0:lmax,0:mmax) :: alm
  !internal
  integer :: nside
  double precision, allocatable :: S(:,:)
  !rmargs :: npix
  !add2py :: npix = len(map0)

  nside = int(dsqrt(npix/12d0))

  allocate(S(0:npix-1,2))
  S(:,1) = map0
  S(:,2) = map1
  call map2alm_spin(nside,lmax,mmax,spin,S,alm)
  deallocate(S)

end subroutine hp_map2alm_spin

subroutine map_mul_lfunc(npix,imap,lmax,lfunc,omap)
!*  Convert map to alm, multiply a function to alm and convert back again to map
!*
!*  Args:
!*    :imap [pix] (double)  : Input map, with bounds (0:12*nside**2-1)
!*    :lfunc [l] (double)   : 1D spectrum to be multiplied to alm, with bounds (0:lmax)
!*
!*  Returns:
!*    :omap [pix] (double)  : Output map, with bounds (0:12*nside**2-1)
!*
  !f2py intent(in) lmax, npix, imap, lfunc
  !f2py intent(out) omap
  !f2py depend(npix) imap, omap
  !f2py depend(lmax) lfunc
  implicit none
  !I/O
  integer, intent(in) :: lmax, npix
  double precision, intent(in), dimension(0:npix-1) :: imap
  double precision, intent(in), dimension(0:lmax) :: lfunc
  double precision, intent(out), dimension(0:npix-1) :: omap
  !internal
  integer :: l, nside
  double complex, allocatable :: alm(:,:,:)
  !rmargs :: lmax, npix
  !add2py :: lmax = len(lfunc) - 1
  !add2py :: npix = len(imap)

  nside = int(dsqrt(npix/12d0))

  allocate(alm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,imap,alm)
  do l = 0, lmax
    alm(1,l,:) = alm(1,l,:) * lfunc(l)
  end do
  call alm2map(nside,lmax,lmax,alm,omap)
  deallocate(alm)

end subroutine map_mul_lfunc

subroutine mulwin(npix,lmax,mmax,alm,win,wlm)
!*  Multiply window to a map obtained from alm
!*
!*  Args:
!*    :alm [l,m] (dcmplx) : Harmonic coefficient to be multiplied at window, with bounds (0:lmax,0:mmax)
!*    :win [pix] (double) : Transformed map, with bounds (0:npix-1)
!*
!*  Returns:
!*    :wlm [l,m] (dcmplx) : Harmonic coefficient of the window-multiplied map, with bounds (0:lmax,0:mmax)
!*
  !f2py intent(in) npix, lmax, mmax, alm, win
  !f2py intent(out) wlm
  !f2py depend(lmax) alm, wlm
  !f2py depend(mmax) alm, wlm
  !f2py depend(npix) win
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, mmax
  double complex, intent(in), dimension(0:lmax,0:mmax) :: alm
  double precision, intent(in), dimension(0:npix-1) :: win
  double complex, intent(out), dimension(0:lmax,0:mmax) :: wlm
  !internal
  integer :: nside
  double precision, allocatable :: map(:)
  double complex, allocatable :: tlm(:,:,:)
  !rmargs :: npix, lmax, mmax
  !add2py :: npix = len(win)
  !add2py :: lmax = len(alm[:,0]) - 1
  !add2py :: mmax = len(alm[0,:]) - 1

  nside = int(dsqrt(npix/12d0))

  allocate(map(0:npix-1),tlm(1,0:lmax,0:mmax))
  tlm(1,:,:) = alm
  call alm2map(nside,lmax,mmax,tlm,map)
  map = win*map
  call map2alm(nside,lmax,mmax,map,tlm)
  wlm = tlm(1,:,:)
  deallocate(map,tlm)

end subroutine mulwin

subroutine mulwin_spin(npix,lmax,mmax,spin,elm,blm,win,wlm)
!*  Multiply window to a map obtained from alm
!*
!*  Args:
!*    :elm [l,m] (dcmplx) : Spin-s E-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
!*    :blm [l,m] (dcmplx) : Spin-s B-like harmonic coefficient to be transformed to a map, with bounds (0:lmax,0:mmax)
!*    :win [pix] (double) : Transformed map, with bounds (0:npix-1)
!*
!*  Args (Optional):
!*    :spin (int)         : Spin of the transform, default to 2
!*
!*  Returns:
!*    :wlm [2,l,m] (dcmplx): Parity-eve/odd harmonic coefficients obtained from the window-multiplied map, with bounds (2,0:lmax,0:mmax)
!*
  !f2py intent(in) npix, lmax, mmax, spin, elm, blm, win
  !f2py intent(out) wlm
  !f2py depend(lmax) elm, blm, wlm
  !f2py depend(mmax) elm, blm, wlm
  !f2py depend(npix) win
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, mmax, spin
  double complex, intent(in), dimension(0:lmax,0:mmax) :: elm, blm
  double precision, intent(in), dimension(0:npix-1) :: win
  double complex, intent(out), dimension(2,0:lmax,0:mmax) :: wlm
  !internal
  integer :: nside
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: alm(:,:,:)
  !rmargs :: npix, lmax, mmax
  !opt4py :: spin = 2
  !add2py :: npix = len(win)
  !add2py :: lmax = len(elm[:,0]) - 1
  !add2py :: mmax = len(elm[0,:]) - 1

  nside = int(dsqrt(npix/12d0))

  allocate(map(0:npix-1,2),alm(2,0:lmax,0:mmax))
  alm(1,:,:) = elm
  alm(2,:,:) = blm
  call alm2map_spin(nside,lmax,mmax,spin,alm,map)
  map(:,1) = win*map(:,1)
  map(:,2) = win*map(:,2)
  call map2alm_spin(nside,lmax,mmax,spin,map,wlm)
  deallocate(map,alm)

end subroutine mulwin_spin

subroutine lm_healpy2healpix(lmpy,almpy,lmax,almpix)
!*  Transform healpy alm to healpix alm
!*
!*  Args:
!*    :lmax (int)           : Maximum multipole of the input/output alm satisfying 2 x lmpy = (lmax+1) x (lmax+2)
!*    :almpy[index] (dcmplx): Healpy alm, with bounds (0:lmpy-1)
!*
!*  Returns:
!*    :almpix [l,m] (dcmplx): Healpix alm, with bounds (0:lmax,0:lmax)
!*
  !f2py intent(in) lmpy, lmax, almpy
  !f2py intent(out) almpix
  !f2py depend(lmpy) almpy
  !f2py depend(lmax) almpix
  implicit none
  integer, intent(in) :: lmpy, lmax
  double complex, intent(in), dimension(0:lmpy-1) :: almpy
  double complex, intent(out), dimension(0:lmax,0:lmax) :: almpix
  integer :: l, m, i
  !rmargs :: lmpy, lmax
  !add2py :: lmpy = len(almpy)
  !add2py :: lmax = int((-3.+numpy.sqrt(1.+8.*lmpy))/2.)

  !i = m*(lmax+1) - m(m-1)/2 + (l-m)
  almpix = 0d0
  if (2*lmpy/=(lmax+1)*(lmax+2)) stop 'error (lm_healpy2healpix): size of almpy is inconsistent with that of almpix'
  do m = 0, lmax
    do l = m, lmax
      i = m*lmax - int(m*(m-1d0)/2d0) + l
      almpix(l,m) = almpy(i)
    end do
  end do

end subroutine lm_healpy2healpix

subroutine cosin_healpix(npix,cosin)
!*  Return cos(theta) as a function of the Healpix pixel index
!*
!*  Args:
!*    :nside (int) : Nside of the desired map
!*
!*  Returns:
!*    :cosin[pix] (double) : cosin(theta), with bounds (0:npix-1)
!*
  !I/O
  !f2py intent(in) npix
  !f2py intent(out) cosin
  !f2py depend(npix) cosin
  implicit none
  integer, intent(in) :: npix
  double precision, intent(out), dimension(0:npix-1) :: cosin
  !internal
  integer :: nside, lmax
  double complex, allocatable :: alm(:,:,:)
  !chargs :: npix -> nside
  !add2py :: npix = 12*nside**2

  nside = int(sqrt(npix/12d0))
  lmax = 2*nside

  allocate(alm(1,0:lmax,0:lmax));  alm=0d0
  alm(1,1,0) = 1d0
  cosin = 0d0
  call alm2map(nside,lmax,lmax,alm,cosin)
  cosin = cosin*dsqrt(pi/0.75d0)
  deallocate(alm)

end subroutine cosin_healpix

subroutine load_defangle_takahashi(fname,npix,theta,phi,verbose)
!*  Read theta and phi coordinates at source plane obtained by Takahashi et al. (2017)
!*
!*  Args:
!*    :fname (str) : file name
!*    :npix (int) : Number of pixels of theta and phi data
!*
!*  Args (optional):
!*    :verbose (bool) : output messages, default to False
!*
!*  Returns:
!*    :theta[pix] (double) : theta, with bounds (0:npix-1)
!*    :phi[pix] (double)   : phi, with bounds (0:npix-1)
!*
  !f2py intent(in) fname, npix, verbose
  !f2py intent(out) theta, phi
  !f2py depend(npix) theta, phi
  implicit none
  character(*), intent(in) :: fname
  integer, intent(in) :: npix
  logical, intent(in) :: verbose
  double precision, intent(out), dimension(0:npix-1) :: theta, phi
  !internal
  integer :: nside, rnpix, n
  !opt4py :: verbose = False

  if (verbose) write(*,*) 'read deflection map'
  
  open(12,file=fname,status='old',form='unformatted')
  read(12) nside, rnpix
  
  if (verbose) write(*,*) nside, rnpix
  if (npix/=rnpix) stop 'input npix is inconsistent with data file'

  if (verbose) write(*,*) 'read data'
  
  read(12) (theta(n),n=0,npix-1)
  read(12) (phi(n),n=0,npix-1)
  
  close(12)

end subroutine load_defangle_takahashi

subroutine polcoord2angle(npix,theta,phi,angle,verbose)
!*  Converting theta and phi coordinates at source plane to deflection angle.
!*  The algorithm is provided by Takashi Hamana and Ryuichi Takahashi.
!*
!*  Args:
!*    :npix (int) : Number of pixels of theta and phi data
!*    :theta[pix] (double) : theta, with bounds (0:npix-1)
!*    :phi[pix] (double)   : phi, with bounds (0:npix-1)
!*
!*  Args (optional):
!*    :verbose (bool) : output messages, default to False
!*
!*  Returns:
!*    :angle[pix,2] (double) : deflection angle vector containing two components, with bounds (0:npix-1,1:2)
!*
  !f2py intent(in) npix, verbose, theta, phi
  !f2py intent(out) angle
  !f2py depend(npix) theta, phi, angle
  implicit none
  integer, intent(in) :: npix
  logical, intent(in) :: verbose
  double precision, intent(in), dimension(0:npix-1) :: theta, phi
  double precision, intent(out), dimension(0:npix-1,2) :: angle
  !internal
  integer :: nside, n
  double precision :: theta0, phi0, deltaphi, cosalp, sinalp, alpha, cosdelta, sindelta
  double precision, allocatable :: theta_i(:), phi_i(:)
  !opt4py :: verbose = False

  if (verbose) write(*,*) 'nside', nside
  if (verbose) write(*,*) 'size:', size(theta), size(phi)
  nside = int(dsqrt(npix/12d0))

  !healpix location
  allocate(theta_i(0:npix-1),phi_i(0:npix-1))
  if (verbose) write(*,*) 'obtain image plane theta/phi'
  do n = 0, npix-1
      call pix2ang_ring(nside,n,theta0,phi0)
      theta_i(n) = theta0
      phi_i(n)   = phi0
  end do

  if (verbose) write(*,*) 'theta/phi to angle'
  do n = 0, npix-1
      deltaphi = phi(n) - phi_i(n)
      cosalp   = dcos(theta_i(n))*dcos(theta(n)) + dsin(theta_i(n))*dsin(theta(n))*dcos(deltaphi)
      if (cosalp.gt.1d0) then
          alpha = 0d0
      else
          alpha = dacos(cosalp)
      end if
      sinalp = dsin(alpha)
      if (sinalp.eq.0d0) then
          angle(n,1) = 0d0
          angle(n,2) = 0d0
      else
          cosdelta = (dcos(theta(n))-dcos(theta_i(n))*cosalp)/(dsin(theta_i(n))*sinalp)
          sindelta = dsin(deltaphi)*dsin(theta(n))/sinalp
          angle(n,1) = -alpha*cosdelta
          angle(n,2) = alpha*sindelta
      end if
  end do
  
  deallocate(theta_i,phi_i)

end subroutine polcoord2angle

subroutine polcoord2angle_alm(nside,lmax,theta,phi,glm,clm,verbose)
!*  Converting theta and phi coordinates at source plane to deflection angle.
!*
!*  Args:
!*    :npix (int) : Number of pixels of theta and phi data
!*    :lmax (int) : Maximum multipole of alms for gradient and curl modes
!*    :theta[pix] (double) : theta, with bounds (0:npix-1)
!*    :phi[pix] (double)   : phi, with bounds (0:npix-1)
!*
!*  Args (optional):
!*    :verbose (bool) : output messages, default to False
!*
!*  Returns:
!*    :glm[l,m] (dcmplx) : gradient mode, with bounds (0:lmax,0:lmax)
!*    :clm[l,m] (dcmplx) : curl mode, with bounds (0:lmax,0:lmax)
!*
  !f2py intent(in) verbose, nside, lmax, theta, phi
  !f2py intent(out) glm, clm
  !f2py depend(lmax) glm, clm
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nside, lmax
  double precision, intent(in), dimension(:) :: theta, phi
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: npix, l
  double precision, allocatable :: angle(:,:)
  double complex, allocatable :: dlm(:,:,:)
  !opt4py :: verbose = False

  npix = 12*nside**2

  if (verbose) then
    if (size(theta)/=npix) stop 'size of theta is not npix'
    if (size(phi)/=npix)   stop 'size of phi is not npix'
  end if

  if (verbose) write(*,*) 'convert to angle'
  allocate(angle(0:npix-1,2))
  call polcoord2angle(npix,theta,phi,angle,verbose)

  if (verbose) write(*,*) 'deflection angle to its alms'
  allocate(dlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,angle,dlm)
  deallocate(angle)

  if (verbose) write(*,*) 'convert to glm and clm'
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
      glm(l,0:l) = dlm(1,l,0:l)/dsqrt(dble(l**2+l))
      clm(l,0:l) = dlm(2,l,0:l)/dsqrt(dble(l**2+l))
  end do
  deallocate(dlm)

end subroutine polcoord2angle_alm

subroutine calc_mfs(bn,nu,lmax,walm,V,nside)
!*  Compute 2D Minkowski functionals from a given alm
!*
!*  Args:
!*    :bn (int)            : Number of nu bins
!*    :nu [bin] (double)   : Nu bins, with bounds (bn)
!*    :lmax (int)          : Maximum multipole of the input walm
!*    :walm [l,m] (dcmplx) : Alm with filtering, possibly divided by the map variance, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional): 
!*    :nside (int): Nside of the intermediate map, default to the closest power of 2 to 3xlmax
!*
!*  Returns:
!*    :V [bin,type] (double): The three Minkowski functionals, V0, V1 and V2, at each nu bin, with bounds (bn,0:2)
!*
  !f2py intent(in) lmax, bn, nu, walm, nside
  !f2py intent(out) V
  !f2py depend(bn) nu, V
  !f2py depend(lmax) walm
  implicit none 
  !I/O
  integer, intent(in) :: lmax, bn
  double precision, intent(in), dimension(bn) :: nu
  double complex, intent(in), dimension(0:lmax,0:lmax) :: walm
  double precision, intent(out), dimension(bn,0:2) :: V
  integer, intent(in) :: nside
  !internal
  integer :: n, i, l, npix, nside_t
  double precision :: dnu, dtt, dtp, dpp, atant, mcV(0:2)
  double precision, allocatable :: cost(:), der0(:), der1(:,:), der2(:,:)
  double complex, allocatable :: alm(:,:,:)
  !opt4py :: nside = 0

  nside_t = nside
  if (nside==0)  nside_t = 2**(int(dlog(dble(3*lmax))/dlog(2d0)))
  npix = 12*nside_t**2
  dnu  = nu(2)-nu(1)

  ! compute derivatives 
  allocate(alm(1,0:lmax,0:lmax),der0(0:npix-1),der1(0:npix-1,2),der2(0:npix-1,3))
  alm(1,:,:) = walm
  call alm2map_der(nside,lmax,lmax,alm,der0,der1,der2)
  deallocate(alm)

  ! cosine
  allocate(cost(0:npix-1))
  call cosin_healpix(npix,cost) 

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

subroutine mock_galaxy_takahashi(fname,zn,ngz,zi,b0,btype,a,b,zm,sz,zbias,gmap)
!*  Compute galaxy overdensity map from dark matter density map of Takahashi et al. (2017). 
!*  The galaxy z distribution is assumed to have the functional form given by Eq.(7) of https://arxiv.org/abs/1810.03346
!*
!*  Args:
!*    :fname (str)         : Filename of density map
!*    :zn (int)            : Number of tomographic bins
!*    :ngz[zn] (double)    : Total number of galaxies at each z-bin
!*    :zi[zn+1] (double)   : Redshift intervals of the tomographic bins
!*
!*  Args(optional): 
!*    :a (double)          : galaxy distribution shape parameter
!*    :b (double)          : galaxy distribution shape parameter
!*    :zm (double)         : mean redshift of galaxy distribution
!*    :b0 (double)         : constant galaxy bias at z=0
!*
!*  Returns:
!*    :gmap [pix,zbin] (double): The galaxy number density map at each zbin
!*
  !f2py intent(in) fname, btype, zn, ngz, zi, b0, a, b, zm, sz, zbias
  !f2py intent(out) gmap
  !f2py depend(zn) ngz, zi, gmap
  implicit none
  character(*), intent(in) :: fname, btype
  integer, intent(in) :: zn
  double precision, intent(in), dimension(1:zn) :: ngz
  double precision, intent(in), dimension(1:zn+1) :: zi
  double precision, intent(in) :: b0, a, b, zm, sz, zbias
  double precision, intent(out), dimension(0:201326591,zn) :: gmap
  !internal
  integer :: npix, nside, kplane, ngal
  integer :: ipix, i, j, k
  real(4) :: zk
  real(4), allocatable :: deltam(:)
  double precision :: z(42), dn(42), omegam, norm, mu, mu_mean
  !opt4py :: a = 2.0
  !opt4py :: b = 1.0
  !opt4py :: zm = 1.0
  !opt4py :: sz = 0.0
  !opt4py :: zbias = 0.0
  !opt4py :: b0 = 1.0
  !opt4py :: btype = 'sqrtz'
 
  ! simulation parameters
  nside  = 4096
  npix   = 12*nside**2
  omegam = 0.279d0

  dn     = 0d0
  gmap   = 0d0

  do j = 1, zn

      write(*,*) 'zi =', j

      ! normalization
      open(12,file=fname,status='old',form='unformatted')
      do k = 1, 42
        read(12) kplane, zk !index of lens plane and redshift
        z(k)  = zk
        if (z(k)<6d0) dn(k) = 5.004d-2*dsqrt(omegam*(1d0+z(k))**3+(1d0-omegam))*nz_SF_scal(z(k),a,b,zm)*pz_SF_scal(z(k),zi(j:j+1),sz,zbias)
      end do
      norm = ngz(j)/sum(dn)

      allocate(deltam(0:npix-1))

      do k = 1, 42

        write(*,*) k
        read(12) kplane
        read(12) deltam(0:npix-1) ! delta^m at the lens plane

        if (z(k)>=6d0) cycle

        mu_mean = norm*dn(k)/dble(npix) ! mean number density of galaxies at kth plane

        do ipix = 0, npix-1
          mu = mu_mean * ( 1d0 + b0*gbias(z(k),btype)*deltam(ipix) )
          if (mu>0.) then
            call poisson(mu,ngal)
          else
            ngal = 0
          end if
          gmap(ipix,j) = gmap(ipix,j) + real(ngal)
        end do

      end do

      close(12)

      deallocate(deltam)

  end do

end subroutine mock_galaxy_takahashi



end module utils


