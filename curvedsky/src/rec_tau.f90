!////////////////////////////////////////////////////!
! Amplitude modulation reconstruction
!////////////////////////////////////////////////////!

module rec_tau
  use alm_tools, only: alm2map, map2alm, alm2map_spin
  use constants, only: iu

  private alm2map, map2alm, alm2map_spin
  private iu

contains 


subroutine qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,alm,nside_t,verbose)
!*  Reconstructing inhomogeneous tau from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output tau alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TT spectrum, with bounds (1:rlmax)
!*    :Tlm1 [l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Tlm2 [l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)      : Nside for the convolution calculation
!*    :verbose (bool)     : Output messages, default to False
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Amplitude modulation alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm1, Tlm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, npix, nside
  double precision, allocatable :: map1(:), map2(:)
  double complex, allocatable :: zlm1(:,:,:), zlm2(:,:,:)
  !opt4py :: nside_t = 0
  !opt4py :: verbose = False

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)   write(*,*) 'calc tau-TT estimator with nside=', nside
  npix = 12*nside**2

  ! alm to map 
  allocate(zlm1(1,0:rlmax,0:rlmax),zlm2(1,0:rlmax,0:rlmax)); zlm1 = 0d0; zlm2=0d0
  do l = rlmin, rlmax
    zlm1(1,l,0:l) = Tlm1(l,0:l)
    zlm2(1,l,0:l) = fC(l)*Tlm2(l,0:l)
  end do 
  allocate(map1(0:npix-1),map2(0:npix-1))
  call alm2map(nside,rlmax,rlmax,zlm1,map1)
  call alm2map(nside,rlmax,rlmax,zlm2,map2)
  deallocate(zlm1,zlm2)

  ! map to alm
  allocate(zlm1(1,0:lmax,0:lmax))
  map1 = map1*map2
  call map2alm(nside,lmax,lmax,map1,zlm1)
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = zlm1(1,l,0:l)
  end do
  deallocate(map1,map2,zlm1)

end subroutine qtt


subroutine qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,alm,nside_t,verbose)
!*  Reconstructing amplitude modulation from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fCE [l] (double)   : EE spectrum, with bounds (0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)      : Nside for the convolution calculation
!*    :verbose (bool)     : Output messages, default to False
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Amplitude modulation alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  !opt4py :: nside_t = 0
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:rlmax) :: fCE
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, nside, npix
  double precision, allocatable :: map(:), A(:,:), A2(:,:)
  double complex, allocatable :: zlm(:,:,:)

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)   write(*,*) 'calc tau-EB estimator with nside=', nside
  npix = 12*nside**2

  ! convolution
  allocate(zlm(2,0:rlmax,0:rlmax),A(0:npix-1,2),A2(0:npix-1,2))
  zlm = 0d0
  do l = rlmin, rlmax
    zlm(2,l,0:l) = Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,zlm,A)
  zlm = 0d0
  do l = rlmin, rlmax
    zlm(1,l,0:l) = fCE(l)*Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,zlm,A2)
  deallocate(zlm)

  allocate(map(0:npix-1))
  map = A(:,1)*A2(:,1)+A(:,2)*A2(:,2)
  deallocate(A,A2)

  allocate(zlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map,zlm)
  deallocate(map)

  ! compute alm
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = zlm(1,l,0:l)
  end do
  deallocate(zlm)

end subroutine qeb


end module rec_tau


