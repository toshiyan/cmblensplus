!////////////////////////////////////////////////////!
! * CMB lensing reconstruction in curvedsky
!////////////////////////////////////////////////////!

module rec_lens
  use alm_tools, only: alm2map, alm2map_spin, map2alm_spin
  use constants, only: iu

  private alm2map, alm2map_spin, map2alm_spin
  private iu

contains 


subroutine qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,glm,clm,nside,gtype)
!*  Reconstructing CMB lensing potential and its curl mode from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TT spectrum, with bounds (0:rlmax)
!*    :Tlm1 [l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Tlm2 [l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int) : Nside for the convolution calculation, default to lmax
!*    :gtype (str) : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential alm, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential) alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  character(1), intent(in), optional :: gtype
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !f2py character(1) :: gtype = ''
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm1, Tlm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, ns, npix
  double precision, allocatable :: at(:), map(:,:), ilk(:)
  double complex, allocatable :: alm1(:,:,:), blm(:,:,:)

  ns = lmax
  if (present(nside)) ns = nside

  allocate(ilk(lmax)); ilk = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  write(*,*) 'calc qTT lens estimator with nside=', ns
  npix = 12*ns**2

  ! compute convolution
  allocate(alm1(1,0:rlmax,0:rlmax))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Tlm1(l,0:l)
  end do 
  allocate(at(0:npix-1))
  call alm2map(nside,rlmax,rlmax,alm1,at)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = fC(l)*Tlm2(l,0:l)*dsqrt(dble((l+1)*l))
  end do 
  allocate(map(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,map)
  map(:,1) = at*map(:,1)
  map(:,2) = at*map(:,2)
  deallocate(at,alm1)

  allocate(blm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,blm)
  deallocate(map)

  ! compute glm and clm
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*blm(1,l,0:l)
    clm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*blm(2,l,0:l)
  end do

  deallocate(ilk,blm)

end subroutine qtt


subroutine qte(lmax,rlmin,rlmax,fC,Tlm,Elm,glm,clm,nside,gtype)
!*  Reconstructing CMB lensing potential and its curl mode from the TE quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TE spectrum, with bounds (0:rlmax)
!*    :Tlm [l,m] (dcmplx) : Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int) : Nside for the convolution calculation, default to lmax
!*    :gtype (str) : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none 
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  character(1), intent(in), optional :: gtype
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !f2py character(1) :: gtype = ''
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in),  dimension(0:rlmax,0:rlmax) :: Tlm, Elm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, ns, npix
  double precision, dimension(:), allocatable :: AT, ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, AE, map
  double complex, dimension(:,:,:), allocatable :: alm1, alm3, blm

  ns = lmax
  if (present(nside)) ns = nside

  allocate(ilk(lmax)); ilk = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  write(*,*) 'calc qTE lens estimator with nside=', ns

  npix = 12*ns**2

  ! convolution
  allocate(alm1(2,0:rlmax,0:rlmax),alm3(2,0:rlmax,0:rlmax),A(0:npix-1,2),A1(0:npix-1,2),A3(0:npix-1,2))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,alm1,A)
  alm1 = 0d0
  alm3 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = fC(l)*Tlm(l,0:l)*dsqrt(dble((l+2)*(l-1)))
    alm3(1,l,0:l) = fC(l)*Tlm(l,0:l)*dsqrt(dble((l-2)*(l+3)))
  end do 
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,A1)
  call alm2map_spin(nside,rlmax,rlmax,3,alm3,A3)
  deallocate(alm1,alm3)

  allocate(alm1(1,0:rlmax,0:rlmax),AT(0:npix-1))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Tlm(l,0:l)
  end do
  call alm2map(nside,rlmax,rlmax,alm1,AT)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax),AE(0:npix-1,2))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Elm(l,0:l)*fC(l)*dsqrt(dble(l*(l+1)))
  end do
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,AE)
  deallocate(alm1)

  allocate(map(0:npix-1,2))
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))  + AT*AE(:,1)*2d0
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1)) + AT*AE(:,2)*2d0
  deallocate(A,A1,A3,AT,AE)

  allocate(blm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,blm)
  deallocate(map)

  ! compute glm and clm
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*blm(1,l,0:l)
    clm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*blm(2,l,0:l)
  end do
  deallocate(blm)

end subroutine qte


subroutine qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,glm,clm,nside,gtype)
!*  Reconstructing CMB lensing potential and its curl mode from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TE spectrum, with bounds (0:rlmax)
!*    :Tlm [l,m] (dcmplx) : Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)  : Nside for the convolution calculation, default to lmax
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none 
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  character(1), intent(in), optional :: gtype
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !f2py character(1) :: gtype = ''
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, ns, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, map
  double complex, dimension(:,:,:), allocatable :: alm1, alm3, zlm

  ns = lmax
  if (present(nside)) ns = nside

  allocate(ilk(lmax)); ilk = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  write(*,*) 'calc qTB lens estimator with nside=', ns
  npix = 12*ns**2

  ! convolution
  allocate(alm1(2,0:rlmax,0:rlmax))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(2,l,0:l) = Blm(l,0:l)
  end do 
  allocate(A(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax),alm3(2,0:rlmax,0:rlmax))
  alm1 = 0d0;  alm3 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = fC(l)*Tlm(l,0:l)*dsqrt(dble((l+2)*(l-1)))
    alm3(1,l,0:l) = fC(l)*Tlm(l,0:l)*dsqrt(dble((l-2)*(l+3)))
  end do 
  allocate(A1(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,A1)
  deallocate(alm1)
  allocate(A3(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npix-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(zlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,zlm)
  deallocate(map)

  ! compute glm and clm
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*zlm(1,l,0:l)
    clm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*zlm(2,l,0:l)
  end do
  deallocate(zlm)

end subroutine qtb


subroutine qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,glm,clm,nside,gtype)
!*  Reconstructing CMB lensing potential and its curl mode from the EE quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : EE spectrum, with bounds (0:rlmax)
!*    :Elm1 [l,m] (dcmplx): 1st inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Elm2 [l,m] (dcmplx): 2nd inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)  : Nside for the convolution calculation, default to lmax
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none 
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  character(1), intent(in), optional :: gtype
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !f2py character(1) :: gtype = ''
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm1, Elm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, ns, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, map
  double complex, dimension(:,:,:), allocatable :: alm, blm

  ns = lmax
  if (present(nside)) ns = nside

  allocate(ilk(lmax)); ilk = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  write(*,*) 'calc qEE lens estimator with nside=', ns
  npix = 12*ns**2

  ! convolution
  allocate(A(0:npix-1,2),alm(2,0:rlmax,0:rlmax))
  alm = 0d0
  do l = rlmin, rlmax
    alm(1,l,0:l) = Elm1(l,0:l)
  end do
  call alm2map_spin(nside,rlmax,rlmax,2,alm,A)
  deallocate(alm)

  allocate(alm(2,0:rlmax,0:rlmax),blm(2,0:rlmax,0:rlmax))
  alm = 0d0;  blm = 0d0
  do l = rlmin, rlmax
    alm(1,l,0:l) = fC(l)*Elm2(l,0:l)*dsqrt(dble((l+2)*(l-1)))
    blm(1,l,0:l) = fC(l)*Elm2(l,0:l)*dsqrt(dble((l-1)*(l+3)))
  end do
  allocate(A1(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,1,alm,A1)
  deallocate(alm)
  allocate(A3(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,3,blm,A3)
  deallocate(blm)

  allocate(map(0:npix-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(alm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,alm)
  deallocate(map)

  ! compute glm and clm 
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l)*dble(l+1))*alm(1,l,0:l)
    clm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l)*dble(l+1))*alm(2,l,0:l)
  end do
  deallocate(alm)

end subroutine qee


subroutine qeb(lmax,rlmin,rlmax,fC,Elm,Blm,glm,clm,nside,gtype)
!*  Reconstructing CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : EE spectrum, with bounds (0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)  : Nside for the convolution calculation, default to lmax
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  character(1), intent(in), optional :: gtype
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !f2py character(1) :: gtype = ''
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, ns, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A,A1,A3,map
  double complex, dimension(:,:,:), allocatable :: alm1,alm3,tlm

  ns = lmax
  if (present(nside)) ns = nside

  allocate(ilk(lmax)); ilk = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  write(*,*) 'calc qEB lens estimator with nside=', ns
  npix = 12*ns**2

  ! convolution
  allocate(alm1(2,0:rlmax,0:rlmax),A(0:npix-1,2))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(2,l,0:l) = Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax),alm3(2,0:rlmax,0:rlmax),A1(0:npix-1,2))
  alm1 = 0d0
  alm3 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = fC(l)*Elm(l,0:l)*dsqrt(dble((l+2)*(l-1)))
    alm3(1,l,0:l) = fC(l)*Elm(l,0:l)*dsqrt(dble((l-2)*(l+3)))
  end do 
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,A1)
  deallocate(alm1)

  allocate(A3(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npix-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(tlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,tlm)
  deallocate(map)

  ! compute glm and clm 
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*tlm(1,l,0:l)
    clm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*tlm(2,l,0:l)
  end do
  deallocate(tlm)

!
!////  if unlensed BB/=0 ////!
!
!  allocate(alm1(0:rlmax,0:rlmax),alm3(0:rlmax,0:rlmax))
!  alm1 = 0d0;  alm3 = 0d0
!  do l = rlmin, rlmax
!    alm1(l,:) = fCBB(l)*Blm(l,:)*llsq(l,-2)
!    alm3(l,:) = fCBB(l)*Blm(l,:)*llsq(l,2)
!  end do 

!  allocate(A(0:npix-1))
!  call elm2map_spin(nside,rlmax,rlmax,2,Elm,A)

!  allocate(A1(0:npix-1))
!  call blm2map_spin(nside,rlmax,rlmax,1,alm1,A1)
!  deallocate(alm1)

!  allocate(A3(0:npix-1))
!  call blm2map_spin(nside,rlmax,rlmax,3,alm3,A3)
!  deallocate(alm3)

!  allocate(map(0:npix-1))
!  map = A*conjg(A1) - conjg(AE)*A3
!  deallocate(A,A1,A3)

!  allocate(tlm(2,0:lmax,0:lmax))
!  call cmplxmap2alm_spin(nside,lmax,lmax,1,map,tlm)
!  deallocate(map)

!  write(*,*) 'compute grad and curl'
!  do l = 1, lmax
!    est(E%g,l,:) = est(E%g,l,:) + 0.5d0*llsq(l,0)*tlm(1,l,:)
!    est(E%c,l,:) = est(E%c,l,:) + 0.5d0*llsq(l,0)*tlm(2,l,:)
!  end do
!  deallocate(tlm)

end subroutine qeb


subroutine qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,glm,clm,nside,gtype)
!*  Reconstructing CMB lensing potential and its curl mode from the BB quadratic estimator
!*
!*  Args:
!*    :lmax (int)          : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipoles of CMB for reconstruction
!*    :fC [l] (double)     : BB spectrum, with bounds (0:rlmax)
!*    :Blm1 [l,m] (dcmplx) : 1st inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm2 [l,m] (dcmplx) : 2nd inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)  : Nside for the convolution calculation, default to lmax
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  integer, intent(in), optional :: nside
  character(1), intent(in), optional :: gtype
  !f2py integer :: nside = lmax
  !docstr :: nside = lmax
  !f2py character(1) :: gtype = ''
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Blm1, Blm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, ns, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, map
  double complex, dimension(:,:,:), allocatable :: alm1, alm3, tlm

  ns = lmax
  if (present(nside)) ns = nside

  allocate(ilk(lmax)); ilk = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  write(*,*) 'calc qBB lens estimator with nside=', ns
  npix = 12*ns**2

  ! convolution
  allocate(alm1(2,0:rlmax,0:rlmax))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(2,l,0:l) = Blm1(l,0:l)
  end do
  allocate(A(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax),alm3(2,0:rlmax,0:rlmax))
  alm1 = 0d0;  alm3 = 0d0
  do l = rlmin, rlmax
    alm1(2,l,0:l) = fC(l)*Blm2(l,0:l)*dsqrt(dble((l+2)*(l-1)))
    alm3(2,l,0:l) = fC(l)*Blm2(l,0:l)*dsqrt(dble((l-2)*(l+3)))
  end do
  allocate(A1(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,A1)
  allocate(A3(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npix-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(tlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,tlm)
  deallocate(map)

  ! compute glm and clm
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*tlm(1,l,0:l)
    clm(l,0:l) = ilk(l)*0.5d0*dsqrt(dble(l*(l+1)))*tlm(2,l,0:l)
  end do
  deallocate(tlm)

end subroutine qbb

end module rec_lens


