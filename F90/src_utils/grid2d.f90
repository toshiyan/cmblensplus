!////////////////////////////////////////////////////!
! * Utils for 2D Grid Analysis
!////////////////////////////////////////////////////!

module grid2d
  use random, only: InitRandom, ranmar, Gaussian1
  use constants,  only: iu, pi, twopi
  use general,  only: filelines, savetxt, loadtxt
  implicit none

  interface spin_weight
    module procedure spin_weight_1darray, spin_weight_2darray
  end interface 

  interface make_lmask
    module procedure make_lmask_int, make_lmask_dble
  end interface 


  private InitRandom, ranmar, gaussian1
  private iu, pi, twopi
  private filelines, savetxt, loadtxt

contains 


!//// Arrays in Fourier plane ////!

function elxy(i,n,D) result(f)
  implicit none
  integer, intent(in) :: i, n
  double precision, intent(in) :: D
  double precision :: f

  f = twopi*dble(i-1-n*0.5d0)/D

end function elxy


function elarray_x(nn,D)  result(f)
! return ell_x in 2D grid
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: f(1:nn(1)*nn(2))

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      f(n) = elxy(i,nn(1),D(1))
      n = n + 1
    end do
  end do

end function elarray_x


function elarray_y(nn,D)  result(f)
! return ell_y in 2D grid
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: f(nn(1)*nn(2))

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      f(n) = elxy(j,nn(2),D(2))
      n = n + 1
    end do
  end do

end function elarray_y


function elarray(nn,D)  result(f)
! return absolute value of multipole in 2D
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: lx, ly, f(nn(1)*nn(2))

  f = 0d0
  n = 1
  do i = 1, nn(1)
    lx = elxy(i,nn(1),D(1))
    do j = 1, nn(2)
      ly = elxy(j,nn(2),D(2))
      if (lx/=0.or.ly/=0)  f(n) = dsqrt(lx**2+ly**2)
      n = n + 1
    end do
  end do

end function elarray


function elarray_inv(nn,D)  result(f)
! * return absolute value of inverse multipole in 2D
! * avoid pixel at ell=0
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: lx, ly, f(nn(1)*nn(2))

  f = 0d0
  n = 1
  do i = 1, nn(1)
    lx = elxy(i,nn(1),D(1))
    do j = 1, nn(2)
      ly = elxy(j,nn(2),D(2))
      if(lx/=0d0.or.ly/=0d0)  f(n) = 1d0/dsqrt(lx**2+ly**2)
      n = n + 1
    end do
  end do

end function elarray_inv


subroutine elarrays_1d(nn,D,elx,ely,els,eli,ei2p)
!* Return elx, ely and absolute value of multipole as 1D array
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(out), optional :: elx(:), ely(:), els(:), eli(:)
  double complex, intent(out), optional :: ei2p(:)
  !internal
  integer :: i, j, n, npix
  double precision :: lx, ly
  double precision, dimension(:), allocatable :: elx0, ely0, els0, eli0
  double complex, allocatable :: ei2p0(:)

  npix = nn(1)*nn(2)
  allocate(elx0(npix),ely0(npix),els0(npix),eli0(npix),ei2p0(npix))
  elx0=0d0;  ely0=0d0;  els0=0d0;  eli0=0d0;  ei2p0=0d0

  n = 0
  do i = 1, nn(1)
    lx = elxy(i,nn(1),D(1))
    do j = 1, nn(2)
      n = n + 1
      ly = elxy(j,nn(2),D(2))
      elx0(n) = lx
      ely0(n) = ly
      if (lx==0d0.and.ly==0d0) cycle
      els0(n)  = dsqrt(lx**2+ly**2)
      eli0(n)  = 1d0/els0(n)
      !compute exp(2i*phi) = cos(2phi) + i*sin(2phi)
      ei2p0(n) = ((lx**2-ly**2)+iu*2d0*lx*ly)/(lx**2+ly**2)
    end do
  end do

  if (present(elx)) elx=elx0
  if (present(ely)) ely=ely0
  if (present(els)) els=els0
  if (present(eli)) eli=eli0
  if (present(ei2p)) ei2p=ei2p0

  deallocate(elx0,ely0,els0,eli0,ei2p0)

end subroutine elarrays_1d


subroutine elarrays_2d(nn,D,elx,ely,els,eli,ei2p)
!* Return elx, ely and absolute value of multipole as 2D array
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(out), optional :: elx(:,:), ely(:,:), els(:,:), eli(:,:)
  double complex, intent(out), optional :: ei2p(:,:)
  !internal
  integer :: i, j
  double precision :: lx, ly
  double precision, allocatable :: elx0(:,:), ely0(:,:), els0(:,:), eli0(:,:)
  double complex, allocatable :: ei2p0(:,:)

!write(*,*) 'E'
  allocate(elx0(nn(1),nn(2)),ely0(nn(1),nn(2)),els0(nn(1),nn(2)),eli0(nn(1),nn(2)),ei2p0(nn(1),nn(2)))
  elx0=0d0;  ely0=0d0;  els0=0d0;  eli0=0d0;  ei2p0=0d0
!write(*,*) 'E'

  do i = 1, nn(1)
    lx = elxy(i,nn(1),D(1))
!write(*,*) 'X'
    do j = 1, nn(2)
      ly = elxy(j,nn(2),D(2))
      elx0(i,j) = lx
      ely0(i,j) = ly
      if (lx==0d0.and.ly==0d0) cycle
      els0(i,j)  = dsqrt(lx**2+ly**2)
      eli0(i,j)  = 1d0/els0(i,j)
      ei2p0(i,j) = ((lx**2-ly**2)+iu*2d0*lx*ly)/(lx**2+ly**2)
    end do
  end do

  if (present(elx)) elx=elx0
  if (present(ely)) ely=ely0
  if (present(els)) els=els0
  if (present(eli)) eli=eli0
  if (present(ei2p)) ei2p=ei2p0

  deallocate(elx0,ely0,els0,eli0,ei2p0)

end subroutine elarrays_2d


subroutine make_lmask_int(nn,D,rL,mask)
  implicit none
  integer, intent(in) :: nn(2), rL(2)
  double precision, intent(in) :: D(2)
  double precision, intent(out) :: mask(:,:)
  ! internal
  integer :: i, j
  double precision, allocatable :: els(:,:), mask_tmp(:,:)

  allocate(els(nn(1),nn(2)),mask_tmp(nn(1),nn(2)))

  call elarrays_2d(nn,D,els=els)

  mask_tmp = 1d0
  do i = 1, nn(1)
    do j = 1, nn(2)
      if (rL(1)>els(i,j).or.els(i,j)>rL(2))  mask_tmp(i,j) = 0d0
    end do
  end do
  mask = mask_tmp
  
  deallocate(els,mask_tmp)

end subroutine make_lmask_int


subroutine make_lmask_dble(nn,D,rL,mask)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), rL(2)
  double precision, intent(out) :: mask(:,:)
  ! internal
  integer :: i, j
  double precision, allocatable :: els(:,:), mask_tmp(:,:)

  allocate(els(nn(1),nn(2)),mask_tmp(nn(1),nn(2)))

  call elarrays_2d(nn,D,els=els)

  mask_tmp = 1d0
  do i = 1, nn(1)
    do j = 1, nn(2)
      if (rL(1)>els(i,j).or.els(i,j)>rL(2))  mask_tmp(i,j) = 0d0
    end do
  end do
  mask = mask_tmp
  
  deallocate(els,mask_tmp)

end subroutine make_lmask_dble


subroutine make_binmask(nn,D,bmin,bmax,bf,bn)
  implicit none
  integer, intent(in) :: nn(2), bn
  double precision, intent(in) :: D(2)
  double precision, intent(in), dimension(bn) :: bmin, bmax
  double precision, intent(out), dimension(bn,nn(1),nn(2)) :: bf
  integer :: b

  bf = 0d0

  do b = 1, bn
    call make_lmask_dble(nn,D,(/bmin(b),bmax(b)/),bf(b,:,:))
  end do

end subroutine make_binmask


subroutine spin_weight_1darray(sw,nn,D,trans,S)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans, S
  double precision, intent(in) :: D(2)
  double complex, intent(out) :: sw(:)
  !internal
  integer :: i, j, n
  double precision :: ss, x, y

  n  = 1
  ss = -trans*S/abs(S)
  do i = 1, nn(1)
    x = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(x/=0d0.or.y/=0d0) sw(n) = sw(n) * ((x**2-y**2)+iu*2d0*ss*x*y)/(x**2+y**2)
      n = n + 1
    end do
  end do

end subroutine spin_weight_1darray


subroutine spin_weight_2darray(sw,nn,D,trans,S)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans, S
  double precision, intent(in) :: D(2)
  double complex, intent(out) :: sw(:,:)
  !internal
  integer :: i, j, n
  double precision :: ss, x, y

  ss = -trans*S/abs(S)
  do i = 1, nn(1)
    x = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(x/=0d0.or.y/=0d0) sw(i,j) = sw(i,j) * ((x**2-y**2)+iu*2d0*ss*x*y)/(x**2+y**2)
    end do
  end do

end subroutine spin_weight_2darray


function WAP(mapsize,s,apfactor,skycut)
!* 1D window function
  implicit none
  !I/O
  double precision, intent(in) :: mapsize, s
  double precision, intent(in), optional :: skycut, apfactor
  !internal
  double precision :: s0, a, Wap, ss, x

  !map size
  a = mapsize*0.5d0
  if(present(skycut)) a = a*skycut

  !apodization length
  s0 = 1d0
  if(present(apfactor)) s0 = apfactor

  !normalized coordinate
  ss = abs(s)/a
  if(s0/=1d0) x = (1d0-ss)/(1d0-s0)
  
  !window function
  if (ss<s0) then
    Wap = 1d0
  else if (ss>=s0.and.ss<1d0) then
    Wap = x - dsin(2*pi*x)/(2*pi) 
  else 
    Wap = 0d0
  end if

end function WAP


end module grid2d


