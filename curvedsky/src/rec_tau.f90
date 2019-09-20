!////////////////////////////////////////////////////!
! Amplitude modulation reconstruction
!////////////////////////////////////////////////////!

module rec_tau
  use alm_tools, only: alm2map, map2alm, alm2map_spin
  use constants, only: iu

  interface qtt
    module procedure qtt_sym, qtt
  end interface

  private alm2map, map2alm, alm2map_spin
  private iu

contains 


subroutine qtt_sym(lmax,rlmin,rlmax,fC,Tlm,alm,nside)
!*  Reconstructing inhomogeneous tau from the temperature quadratic estimator, assuming Tlm1=Tlm2
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output tau alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TT spectrum, with bounds (1:rlmax)
!*    :Tlm [l,m] (dcmplx) : Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)        : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Amplitude modulation alm, with bounds (0:lmax,0:lmax)
!*

  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, npix, ns
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: zlm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  write(*,*) 'calc qTT estimator with nside=', ns
  npix = 12*ns**2

  ! alm to map 
  allocate(zlm(2,0:rlmax,0:rlmax)); alm = 0d0
  do l = rlmin, rlmax
    zlm(1,l,0:l) = Tlm(l,0:l)
    zlm(2,l,0:l) = fC(l)*Tlm(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rlmax,rlmax,zlm(1:1,:,:),map(1,:))
  call alm2map(nside,rlmax,rlmax,zlm(2:2,:,:),map(2,:))
  deallocate(zlm)

  ! map to alm
  allocate(zlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map(1,:)*map(2,:),zlm)
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = zlm(1,l,0:l)
  end do
  deallocate(map,zlm)

end subroutine qtt_sym


subroutine qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,alm,nside)
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
!*    :nside (int)        : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Amplitude modulation alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm1, Tlm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, npix, ns
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: zlm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  write(*,*) 'calc qTT estimator with nside=', ns
  npix = 12*ns**2

  ! alm to map 
  allocate(zlm(2,0:rlmax,0:rlmax)); alm = 0d0
  do l = rlmin, rlmax
    zlm(1,l,0:l) = Tlm1(l,0:l)
    zlm(2,l,0:l) = fC(l)*Tlm2(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rlmax,rlmax,zlm(1:1,:,:),map(1,:))
  call alm2map(nside,rlmax,rlmax,zlm(2:2,:,:),map(2,:))
  deallocate(zlm)

  ! map to alm
  allocate(zlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map(1,:)*map(2,:),zlm)
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = zlm(1,l,0:l)
  end do
  deallocate(map,zlm)

end subroutine qtt


subroutine qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,alm,nside)
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
!*    :nside (int)        : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Amplitude modulation alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  double precision, intent(in), dimension(0:rlmax) :: fCE
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, ns, npix
  double precision, allocatable :: map(:), A(:,:), A2(:,:)
  double complex, allocatable :: zlm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  write(*,*) 'calc qEB estimator with nside=', ns
  npix = 12*ns**2

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


subroutine oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,alm,nside)
!*  Reconstructing amplitude modulation by the odd EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fEB [l] (double)   : EB spectrum, with bounds (0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)  : Nside for the convolution calculation, default to lmax
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Reconstructed alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  double precision, intent(in), dimension(0:rlmax) :: fEB
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, ns, npix
  double precision, allocatable :: map(:), A(:,:), A2(:,:)
  double complex, allocatable :: xlm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  write(*,*) 'calc qEB rotation estimator with nside=', ns
  npix = 12*ns**2

  ! convolution 1
  allocate(xlm(2,0:rlmax,0:rlmax),A(0:npix-1,2),A2(0:npix-1,2))
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(1,l,0:l) = Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A)
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(2,l,0:l) = -fEB(l)*Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A2)

  deallocate(xlm)

  allocate(map(0:npix-1))
  map = A(:,1)*A2(:,2)-A(:,2)*A2(:,1)
  deallocate(A,A2)

  ! convolution 2
  allocate(xlm(2,0:rlmax,0:rlmax),A(0:npix-1,2),A2(0:npix-1,2))
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(1,l,0:l) = fEB(l)*Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A2)
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(2,l,0:l) = Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A)
  deallocate(xlm)

  map = map + A(:,1)*A2(:,2)-A(:,2)*A2(:,1)
  deallocate(A,A2)

  allocate(xlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map,xlm)
  deallocate(map)

  ! compute alm 
  write(*,*) 'compute reconstructed fields'
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = xlm(1,l,0:l)
  end do
  deallocate(xlm)

end subroutine oeb


end module rec_tau


