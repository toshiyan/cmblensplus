!////////////////////////////////////////////////////!
! Inhomogeneous tau reconstruction
!////////////////////////////////////////////////////!

module rec_tau
  use alm_tools, only: alm2map, map2alm
  use constants, only: iu

  interface qtt
    module procedure qtt_sym, qtt
  end interface

  private alm2map, map2alm
  private iu

contains 


subroutine qtt_sym(lmax,rlmin,rlmax,fC,Tlm,xlm,nside)
!*  Reconstructing inhomogeneous tau from the temperature quadratic estimator, assuming Tlm1=Tlm2
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output tau alms
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : Temperature angular power spectrum, with bounds (1:rlmax)
!*    :Tlm [l,m] (dcmplx) : Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)       : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :xlm [l,m] (dcmplx) : Inhomogeneous tau alm, with bounds (0:lmax,0:lmax)
!*

  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
!f2py integer :: nside = lmax
  double precision, intent(in), dimension(rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: xlm
  !internal
  integer :: l, npix, ns
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: alm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  write(*,*) 'calc qTT tau estimator, nside = ', ns
  npix = 12*ns**2

  ! alm to map 
  allocate(alm(2,0:rlmax,0:rlmax)); alm = 0d0
  do l = rlmin, rlmax
    alm(1,l,0:l) = Tlm(l,0:l)
    alm(2,l,0:l) = fC(l)*Tlm(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rlmax,rlmax,alm(1:1,:,:),map(1,:))
  call alm2map(nside,rlmax,rlmax,alm(2:2,:,:),map(2,:))
  deallocate(alm)

  ! map to alm
  allocate(alm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map(1,:)*map(2,:),alm)
  xlm = 0d0
  do l = 1, lmax
    xlm(l,0:l) = alm(1,l,0:l)
  end do
  deallocate(map,alm)

end subroutine qtt_sym


subroutine qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,xlm,nside)
!*  Reconstructing inhomogeneous tau from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output tau alms
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : Temperature angular power spectrum, with bounds (1:rlmax)
!*    :Tlm1 [l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Tlm2 [l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)       : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :xlm [l,m] (dcmplx) : Inhomogeneous tau alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
!f2py integer :: nside = lmax
  double precision, intent(in), dimension(rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm1, Tlm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: xlm
  !internal
  integer :: l, npix, ns
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: alm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  write(*,*) 'calc qTT tau estimator'
  npix = 12*ns**2

  ! alm to map 
  allocate(alm(2,0:rlmax,0:rlmax)); alm = 0d0
  do l = rlmin, rlmax
    alm(1,l,0:l) = Tlm1(l,0:l)
    alm(2,l,0:l) = fC(l)*Tlm2(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rlmax,rlmax,alm(1:1,:,:),map(1,:))
  call alm2map(nside,rlmax,rlmax,alm(2:2,:,:),map(2,:))
  deallocate(alm)

  ! map to alm
  allocate(alm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map(1,:)*map(2,:),alm)
  xlm = 0d0
  do l = 1, lmax
    xlm(l,0:l) = alm(1,l,0:l)
  end do
  deallocate(map,alm)

end subroutine qtt


end module rec_tau


