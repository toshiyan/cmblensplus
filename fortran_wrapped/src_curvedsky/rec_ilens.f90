!////////////////////////////////////////////////////!
! * Imaginary CMB lensing reconstruction in curvedsky
!////////////////////////////////////////////////////!

! need update

module rec_ilens
  use alm_tools, only: alm2map, alm2map_spin, map2alm_spin
  use constants, only: iu

  private alm2map, alm2map_spin, map2alm_spin
  private iu

contains 


subroutine qte(lmax,rlmin,rlmax,fC,Tlm,Elm,glm,clm,nside_t,gtype,verbose)
!*  Reconstructing imaginary CMB lensing potential and its curl mode from the TE quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TB spectrum, with bounds (0:rlmax)
!*    :Tlm [l,m] (dcmplx) : Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int) : Nside for the convolution calculation
!*    :gtype (str) : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*    :verbose (bool) : Output messages, default to False
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none 
  !f2py intent(in) verbose, lmax, rlmin, rlmax, nside_t, gtype, fC, Tlm, Elm
  !f2py intent(out) glm, clm
  !f2py depend(rlmax) fC, Tlm, Elm
  !f2py depend(lmax) glm, clm
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  character(1), intent(in) :: gtype
  !opt4py :: nside_t = 0
  !opt4py :: gtype = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in),  dimension(0:rlmax,0:rlmax) :: Tlm, Elm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, nside, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, map
  double complex, dimension(:,:,:), allocatable :: alm1, alm3, blm

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc qTE ilens estimator with nside=', nside
  npix = 12*nside**2

  allocate(ilk(lmax)); ilk = 1d0
  if (gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

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

  allocate(map(0:npix-1,2))
  map(:,1) = -A(:,1)*(A1(:,2)-A3(:,2)) + A(:,2)*(A1(:,1)-A3(:,1))
  map(:,2) = -A(:,1)*(A1(:,1)+A3(:,1)) - A(:,2)*(A1(:,2)+A3(:,2))
  deallocate(A,A1,A3)

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
  deallocate(blm)

end subroutine qte

subroutine qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,glm,clm,nside_t,gtype,verbose)
!*  Reconstructing imaginary CMB lensing potential and its curl mode from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : TE spectrum, with bounds (0:rlmax)
!*    :Tlm [l,m] (dcmplx) : Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)  : Nside for the convolution calculation
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*    :verbose (bool) : Output messages, default to False
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none 
  !f2py intent(in) verbose, lmax, rlmin, rlmax, nside_t, gtype, fC, Tlm, Blm
  !f2py intent(out) glm, clm
  !f2py depend(rlmax) fC, Tlm, Blm
  !f2py depend(lmax) glm, clm
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  character(1), intent(in) :: gtype
  !opt4py :: nside_t = 0
  !opt4py :: gtype = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, nside, npix
  double precision, dimension(:), allocatable :: ilk, AT
  double precision, dimension(:,:), allocatable :: A, A1, A3, map, AB
  double complex, dimension(:,:,:), allocatable :: alm1, alm3, zlm

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc qTB ilens estimator with nside=', nside
  npix = 12*nside**2

  allocate(ilk(lmax)); ilk = 1d0
  if (gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

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

  allocate(alm1(1,0:rlmax,0:rlmax),AT(0:npix-1))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Tlm(l,0:l)
  end do
  call alm2map(nside,rlmax,rlmax,alm1,AT)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax),AB(0:npix-1,2))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(2,l,0:l) = Blm(l,0:l)*fC(l)*dsqrt(dble(l*(l+1)))
  end do
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,AB)
  deallocate(alm1)

  allocate(map(0:npix-1,2))
  !map = A*conjg(A1) + conjg(A)*A3
  map(:,1) = -A(:,1)*(A1(:,2)-A3(:,2)) + A(:,2)*(A1(:,1)-A3(:,1)) + 2d0*AT*AB(:,2)
  map(:,2) = -A(:,1)*(A1(:,1)+A3(:,1)) - A(:,2)*(A1(:,2)+A3(:,2)) - 2d0*AT*AB(:,1)
  deallocate(A,A1,A3,AT,AB)

  allocate(zlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,zlm)
  deallocate(map)

  ! compute glm and clm
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*zlm(1,l,0:l)
    clm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*zlm(2,l,0:l)
  end do
  deallocate(zlm)

end subroutine qtb

subroutine qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,glm,clm,nside_t,gtype,verbose)
!*  Reconstructing imaginary CMB lensing potential and its curl mode from the EE quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : EE spectrum, with bounds (0:rlmax)
!*    :Elm1 [l,m] (dcmplx): 1st inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Elm2 [l,m] (dcmplx): 2nd inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)  : Nside for the convolution calculation
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*    :verbose (bool) : Output messages, default to False
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none 
  !f2py intent(in) verbose, lmax, rlmin, rlmax, nside_t, gtype, fC, Elm1, Elm2
  !f2py intent(out) glm, clm
  !f2py depend(rlmax) fC, Elm1, Elm2
  !f2py depend(lmax) glm, clm
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  character(1), intent(in) :: gtype
  !opt4py :: nside_t = 0
  !opt4py :: gtype = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm1, Elm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, nside, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, map
  double complex, dimension(:,:,:), allocatable :: alm, blm

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc qEE ilens estimator with nside=', nside
  npix = 12*nside**2

  allocate(ilk(lmax)); ilk = 1d0
  if (gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

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
  map(:,1) = -A(:,1)*(A1(:,2)-A3(:,2)) + A(:,2)*(A1(:,1)-A3(:,1))
  map(:,2) = -A(:,1)*(A1(:,1)+A3(:,1)) - A(:,2)*(A1(:,2)+A3(:,2))
  deallocate(A,A1,A3)

  allocate(alm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,alm)
  deallocate(map)

  ! compute glm and clm 
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*dsqrt(dble(l)*dble(l+1))*alm(1,l,0:l)
    clm(l,0:l) = ilk(l)*dsqrt(dble(l)*dble(l+1))*alm(2,l,0:l)
  end do
  deallocate(alm)

end subroutine qee

subroutine qeb(lmax,rlmin,rlmax,fC,Elm,Blm,glm,clm,nside_t,gtype,verbose)
!*  Reconstructing imaginary CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)    : EE spectrum, with bounds (0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)  : Nside for the convolution calculation
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*    :verbose (bool) : Output messages, default to False
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none
  !f2py intent(in) verbose, lmax, rlmin, rlmax, nside_t, gtype, fC, Elm, Blm
  !f2py intent(out) glm, clm
  !f2py depend(rlmax) fC, Elm, Blm
  !f2py depend(lmax) glm, clm
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  character(1), intent(in) :: gtype
  !opt4py :: nside_t = 0
  !opt4py :: gtype = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, nside, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A,A1,A3,map
  double complex, dimension(:,:,:), allocatable :: alm1,alm3,tlm

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc qEB ilens estimator with nside=', nside
  npix = 12*nside**2

  allocate(ilk(lmax)); ilk = 1d0
  if (gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

  ! convolution
  ! 1st part
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
  map(:,1) = -A(:,1)*(A1(:,2)-A3(:,2)) + A(:,2)*(A1(:,1)-A3(:,1))
  map(:,2) = -A(:,1)*(A1(:,1)+A3(:,1)) - A(:,2)*(A1(:,2)+A3(:,2))
  deallocate(A,A1,A3)

  ! 2nd part
  allocate(alm1(2,0:rlmax,0:rlmax),A(0:npix-1,2))
  alm1 = 0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:rlmax,0:rlmax),alm3(2,0:rlmax,0:rlmax),A1(0:npix-1,2))
  alm1 = 0d0
  alm3 = 0d0
  do l = rlmin, rlmax
    alm1(2,l,0:l) = fC(l)*Blm(l,0:l)*dsqrt(dble((l+2)*(l-1)))
    alm3(2,l,0:l) = fC(l)*Blm(l,0:l)*dsqrt(dble((l-2)*(l+3)))
  end do 
  call alm2map_spin(nside,rlmax,rlmax,1,alm1,A1)
  deallocate(alm1)

  allocate(A3(0:npix-1,2))
  call alm2map_spin(nside,rlmax,rlmax,3,alm3,A3)
  deallocate(alm3)

  map(:,1) = map(:,1) + A(:,1)*(A1(:,2)-A3(:,2)) - A(:,2)*(A1(:,1)-A3(:,1))
  map(:,2) = map(:,2) + A(:,1)*(A1(:,1)+A3(:,1)) + A(:,2)*(A1(:,2)+A3(:,2))
  deallocate(A,A1,A3)

  allocate(tlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,tlm)
  deallocate(map)

  ! compute glm and clm 
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*tlm(1,l,0:l)
    clm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*tlm(2,l,0:l)
  end do
  deallocate(tlm)

end subroutine qeb

subroutine qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,glm,clm,nside_t,gtype,verbose)
!*  Reconstructing imaginary CMB lensing potential and its curl mode from the BB quadratic estimator
!*
!*  Args:
!*    :lmax (int)          : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipoles of CMB for reconstruction
!*    :fC [l] (double)     : BB spectrum, with bounds (0:rlmax)
!*    :Blm1 [l,m] (dcmplx) : 1st inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm2 [l,m] (dcmplx) : 2nd inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)  : Nside for the convolution calculation
!*    :gtype (str)  : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*    :verbose (bool) : Output messages, default to False
!*
!*  Returns:
!*    :glm [l,m] (dcmplx) : CMB lensing potential, with bounds (0:lmax,0:lmax)
!*    :clm [l,m] (dcmplx) : Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)
!*
  implicit none
  !f2py intent(in) verbose, lmax, rlmin, rlmax, nside_t, gtype, fC, Blm1, Blm2
  !f2py intent(out) glm, clm
  !f2py depend(rlmax) fC, Blm1, Blm2
  !f2py depend(lmax) glm, clm
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  character(1), intent(in) :: gtype
  !opt4py :: nside_t = 0
  !opt4py :: gtype = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:rlmax) :: fC
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Blm1, Blm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: glm, clm
  !internal
  integer :: l, nside, npix
  double precision, dimension(:), allocatable :: ilk
  double precision, dimension(:,:), allocatable :: A, A1, A3, map
  double complex, dimension(:,:,:), allocatable :: alm1, alm3, tlm

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc qBB ilens estimator with nside=', nside
  npix = 12*nside**2

  allocate(ilk(lmax)); ilk = 1d0
  if (gtype=='k') then
    do l = 1, lmax
      ilk(l) = 2d0/dble(l*(l+1))
    end do
  end if

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
  map(:,1) = A(:,1)*(A1(:,2)-A3(:,2)) - A(:,2)*(A1(:,1)-A3(:,1))
  map(:,2) = A(:,1)*(A1(:,1)+A3(:,1)) + A(:,2)*(A1(:,2)+A3(:,2))
  deallocate(A,A1,A3)

  allocate(tlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,1,map,tlm)
  deallocate(map)

  ! compute glm and clm
  glm = 0d0
  clm = 0d0
  do l = 1, lmax
    glm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*tlm(1,l,0:l)
    clm(l,0:l) = ilk(l)*dsqrt(dble(l*(l+1)))*tlm(2,l,0:l)
  end do
  deallocate(tlm)

end subroutine qbb


end module rec_ilens



