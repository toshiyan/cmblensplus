!////////////////////////////////////////////////////!
! * Delensing in Fullsky
!////////////////////////////////////////////////////!

module delens
  use alm_tools, only: alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  use constants, only: iu

  private alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  private iu

contains 


subroutine lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,lBlm,nside_t,gtype)
!*  Computing lensing B mode as a convolution of wiener-filtered E-mode and lensing potential
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing B-mode alm
!*    :elmin (int)        : Minimum multipole of wiener-filtered E-mode alm
!*    :elmax (int)        : Maximum multipole of wiener-filtered E-mode alm
!*    :plmin (int)        : Minimum multipole of wiener-filtered lensing potential alm
!*    :plmax (int)        : Maximum multipole of wiener-filtered lensing potential alm
!*    :wElm [l,m] (dcmplx): Wiener-filtered E-mode alm, with bounds (0:elmax,0:elmax)
!*    :wplm [l,m] (dcmplx): Wiener-filtered lensing potential (or kappa) alm, with bounds (0:plmax,0:plmax)
!*
!*  Args(optional):
!*    :nside_t (int)      : Nside for the convolution calculation
!*    :gtype (str)        : Type of input wplm ('p'=phi or 'k'=kappa), default to 'p' (phi). 
!*
!*  Returns:
!*    :lBlm [l,m] (dcmplx): Lensing B-mode alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !f2py intent(in) lmax, elmin, elmax, plmin, plmax, wElm, wplm!wienerfilteredphialm, nside_t, gtype
  !f2py intent(out) lBlm
  !f2py depend(elmax) wElm
  !f2py depend(plmax) wplm!wienerfilteredphialm
  !f2py depend(lmax) lBlm
  !I/O
  integer, intent(in) :: lmax, elmin, elmax, plmin, plmax
  double complex, intent(in), dimension(0:elmax,0:elmax)  :: wElm
  double complex, intent(in), dimension(0:plmax,0:plmax)  :: wplm !wiener filtered phi alm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: lBlm
  !optional
  integer, intent(in) :: nside_t
  character(1), intent(in) :: gtype
  !opt4py :: nside_t = 0
  !opt4py :: gtype = 'p'
  !internal
  integer :: l, m, npix, nside
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A1,A3,A,map
  double complex, dimension(:,:,:), allocatable :: alm

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(max(elmax,plmax)))/dlog(2d0)))
  npix = 12*nside**2

  allocate(ilk(plmax)); ilk = 1d0
  if (gtype=='k') then
    do l = 1, plmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  allocate(A1(0:npix-1,2),A3(0:npix-1,2),A(0:npix-1,2))

  allocate(alm(2,0:elmax,0:elmax)); alm = 0d0
  do l = elmin, elmax
    alm(1,l,:) = WElm(l,:)*dsqrt(dble((l+2)*(l-1))*0.5)
  end do 
  call alm2map_spin(nside,elmax,elmax,1,alm,A1)

  alm = 0d0
  do l = elmin, elmax
    alm(1,l,:) = WElm(l,:)*dsqrt(dble((l-2)*(l+3))*0.5)
  end do 
  call alm2map_spin(nside,elmax,elmax,3,alm,A3)
  deallocate(alm)

  allocate(alm(2,0:plmax,0:plmax)); alm=0d0
  alm = 0d0
  do l = plmin, plmax
    alm(1,l,:) = wplm(l,:)*dsqrt(dble(l*(l+1))*0.5)*ilk(l)
  end do 
  call alm2map_spin(nside,plmax,plmax,1,alm,A)
  deallocate(alm)

  !convolution
  allocate(map(0:npix-1,2))
  !map = A1*A - A3*conjg(A)
  !    = (rA1+iu*iA1)*(rA+iu*iA) - (rA3+iu*iA3)*(rA-iu*iA)
  !    = rA*(rA1-rA3+iu*(iA1-iA3)) + iA*(iu*rA1-iA1+iu*rA3-iA3)
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) - A(:,2)*(A1(:,2)+A3(:,2))
  map(:,2) = A(:,1)*(A1(:,2)-A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A1,A3,A)

  allocate(alm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,2,map,alm)
  deallocate(map)

  lBlm = alm(2,:,:)
  deallocate(alm)

end subroutine lensingb

subroutine shiftvec(npix,lmax,plm,beta,nremap)
!*  Return the anti deflection vector, beta, at the Healpix pixel for the delensing where 
!*
!*    beta(n) + alphaiw(n+beta(n)) = 0
!*
!*  and alphaw is the filtered lensing deflection vector (see arXiv:1701.01712).
!*
!*  Args:
!*    :nside (int)          : Nside of output shift vector
!*    :lmax (int)           : Maximum multipole of the input plm
!*    :plm [l,m] (dcmplx)   : Wiener-filtered lensing potential alm, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :nremap (int)         : Number of iteration for computing the shift vector
!*
!*  Returns:
!*    :beta [pix,2] (double): 2D shift vector, with bounds (0:npix-1,1:2)
!*
  implicit none
  !f2py intent(in) npix, lmax, plm, nremap
  !f2py intent(out) beta
  !f2py depend(lmax) plm
  !f2py depend(npix) beta
  !I/O
  integer, intent(in) :: npix, lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: plm
  double precision, intent(out), dimension(0:npix-1,2) :: beta
  !optional
  integer, intent(in) :: nremap
  !optfpy :: nremap = 3
  !internal
  integer :: nside, j
  double precision, allocatable :: map(:), alpha(:,:), dalpha(:,:)
  double complex, allocatable :: plm0(:,:,:)
  !replace
  !chargs :: npix -> nside
  !add2py :: npix = 12*nside**2
 
  nside = int(dsqrt(npix/12d0))

  allocate(plm0(1,0:lmax,0:lmax),map(0:npix-1),alpha(0:npix-1,2),dalpha(0:npix-1,3))

  plm0(1,:,:) = plm(:,:)
  call alm2map_der(nside,lmax,lmax,plm0,map,alpha,dalpha)

  !1st order itereation
  beta = 0d0
  do j = 1, nremap
    beta(:,1) = alpha(:,1) - dalpha(:,1)*beta(:,1) - dalpha(:,2)*beta(:,2)
    beta(:,2) = alpha(:,2) - dalpha(:,2)*beta(:,1) - dalpha(:,3)*beta(:,2)
  end do
  beta = - beta

  deallocate(plm0,map,alpha,dalpha)

end subroutine shiftvec

subroutine phi2grad(npix,lmax,plm,grad)
!*  Return the deflection vector, grad, at the Healpix pixel
!*
!*  Args:
!*    :nside (int)          : Nside of output deflection vector
!*    :lmax (int)           : Maximum multipole of the input plm/clm
!*    :plm [l,m] (dcmplx)   : Lensing potential alm, with bounds (0:lmax,0:lmax)
!*
!*  Returns:
!*    :grad [pix,2] (double): 2D deflection vector, with bounds (0:npix-1,1:2)
!*
  implicit none
  !f2py intent(in) npix, lmax, plm
  !f2py intent(out) grad
  !f2py depend(lmax) plm
  !f2py depend(npix) grad
  !I/O
  integer, intent(in) :: npix, lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: plm
  double precision, intent(out), dimension(0:npix-1,2) :: grad
  !internal
  integer :: nside
  double precision, allocatable :: map(:)
  double complex, allocatable :: plm0(:,:,:)
  !replace
  !chargs :: npix -> nside
  !add2py :: npix = 12*nside**2

  nside = int(dsqrt(npix/12d0))

  allocate(plm0(1,0:lmax,0:lmax),map(0:npix-1))
  plm0(1,:,:) = plm(:,:)
  call alm2map_der(nside,lmax,lmax,plm0,map,grad)
  deallocate(plm0,map)

end subroutine phi2grad



end module delens



