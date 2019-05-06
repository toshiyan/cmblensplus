!////////////////////////////////////////////////////!
! * Delensing in Fullsky
!////////////////////////////////////////////////////!

module delens
  use alm_tools, only: alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  use constants, only: iu
  use remap_cmb, only: simple_remapping

  private alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  private iu
  private simple_remapping

contains 


subroutine lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,lBlm,nside)
!*  Computing lensing B mode as a convolution of wiener-filtered E-mode and lensing potential
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing B-mode alm
!*    :elmin (int)        : Minimum multipole of wiener-filtered E-mode alm
!*    :elmax (int)        : Maximum multipole of wiener-filtered E-mode alm
!*    :plmin (int)        : Minimum multipole of wiener-filtered lensing potential alm
!*    :plmax (int)        : Maximum multipole of wiener-filtered lensing potential alm
!*    :wElm [l,m] (dcmplx): Wiener-filtered E-mode alm, with bounds (0:elmax,0:elmax)
!*    :wplm [l,m] (dcmplx): Wiener-filtered lensing potential alm, with bounds (0:plmax,0:plmax)
!*
!*  Args(optional):
!*    :nside (int)        : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :lBlm [l,m] (dcmplx): Lensing B-mode alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, elmin, elmax, plmin, plmax
  double complex, intent(in), dimension(0:elmax,0:elmax)  :: wElm
  double complex, intent(in), dimension(0:plmax,0:plmax)  :: wplm !wiener filtered phi alm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: lBlm
  integer, intent(in), optional :: nside
!f2py integer :: nside = lmax
  !internal
  integer :: l, m, npix, ns
  double precision, dimension(:,:), allocatable :: A1,A3,A,map
  double complex, dimension(:,:,:), allocatable :: alm

  ns = lmax
  if (present(nside)) ns = nside
  npix = 12*ns**2

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
    alm(1,l,:) = wplm(l,:)*dsqrt(dble(l*(l+1))*0.5)
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
!*    :npix (int)           : Pixel number of output shift vector
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
  !I/O
  integer, intent(in) :: npix, lmax
  integer, intent(in), optional :: nremap
  double complex, intent(in), dimension(0:lmax,0:lmax) :: plm
  double precision, intent(out), dimension(0:npix-1,2) :: beta
!f2py integer :: nremap = 3
  !internal
  integer :: nside, j
  double precision, allocatable :: map(:), alpha(:,:), dalpha(:,:)
  double complex, allocatable :: plm0(:,:,:)

  !npix = 12*nside**2
  nside = int(dsqrt(npix/12d0))
  write(*,*) 'Nside =', nside

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


subroutine remap_tp(npix,lmax,beta,alm_in,alm_re)
!*  Remapping CMB temperaure and polarization with a given shift vector, grad, based on a simple implementation of LensPix. 
!*  alm[0,:,:] is temperature, alm[1,:,:] is E mode, and alm[2,:,:] is B mode.
!*
!*  Args:
!*    :npix (int)                : Pixel number of output shift vector
!*    :lmax (int)                : Maximum multipole of the input plm
!*    :beta [pix,2] (double)     : 2D shift vector, with bounds (0:npix-1,1:2)
!*    :alm_in [TEB,l,m] (dcmplx) : Input T/E/B alms to be remapped, with bounds (1:3,0:lmax,0:lmax).
!*
!*  Returns:
!*    :alm_re [TEB,l,m] (dcmplx) : Remapped T/E/B alms, with bounds (1:3,0:lmax,0:lmax). 
!*
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax
  double precision, intent(in), dimension(0:npix-1,2) :: beta
  double complex, intent(in), dimension(3,0:lmax,0:lmax) :: alm_in
  double complex, intent(out), dimension(3,0:lmax,0:lmax) :: alm_re
  !internal
  integer :: nside
  real, allocatable :: tqu(:,:)
  complex, allocatable :: alm0(:,:,:), dvec(:)
  double precision, allocatable :: S(:,:)

  !npix = 12*nside**2
  nside = int(dsqrt(npix/12d0))
  write(*,*) 'Nside =', nside

  !if (12*nside**2/=npix) stop 'error (remaptp): inconsistency in nside and npix'

  allocate(alm0(3,0:lmax,0:lmax),dvec(0:npix-1),tqu(0:npix-1,3))
  alm0 = alm_in !double to single
  dvec = cmplx(beta(1:npix,1),beta(1:npix,2))
  call simple_remapping(nside,lmax,alm0,dvec,tqu,'')
  deallocate(alm0,dvec)

  allocate(S(0:npix-1,3))
  S = tqu !single to double
  call map2alm(nside,lmax,lmax,S(:,1),alm_re(1:1,:,:))
  call map2alm_spin(nside,lmax,lmax,2,S(:,2:3),alm_re(2:3,:,:))
  deallocate(S)

end subroutine remap_tp


end module delens


