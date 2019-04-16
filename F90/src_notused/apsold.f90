!////////////////////////////////////////////////////!
! * Angular power spectrum calculation
!////////////////////////////////////////////////////!

module pstool
  use constants, only: pi, TT, TE, EE, BB, dd, Td, Ed, oo
  use general,   only: splint, spline, FileColumns, FileLines, savetxt, linspace, interp_lin, check_positive
  implicit none

  private pi, TT, TE, EE, BB, dd, Td, Ed, oo
  private splint, spline, FileColumns, FileLines, savetxt, linspace, interp_lin, check_positive

contains

#ifdef all 

subroutine angle_average(nn,alm,Al,nmax,label,el)
  implicit none
  !subroutine for computing power spectrum from Fourier grids
  !I/O
  integer, intent(in) :: nn(2), nmax, label(:)
  double precision, intent(in) :: el(:)
  double precision, intent(out) :: Al(:)
  double complex, intent(in) :: alm(:)
  !internal
  integer :: i, m, n
  double precision :: num(nmax)

  !generate [ l(i), Al(i), Cl(i) ] table 
  num = 0
  Al = 0.d0
  do i = 1, nmax 
    do n = 1, nn(1)*nn(2)
      if(label(n)==i) then
        num(i) = num(i) + 1
        Al(i) = Al(i) + alm(n) 
      end if
      Al(i) = Al(i)/dble(num(i))
    end do
  end do

end subroutine angle_average

subroutine power_labeling(el,els,nmax,label)
!* find independent els in 2D grid
  implicit none 
  !I/O
  integer, intent(in) :: el(2)
  double precision, intent(in) :: els(:)
  integer, intent(out) :: nmax, label(:)
  !internal
  integer :: npix, m, n, indep

  npix = size(els)
  nmax = 0

  do n = 1, npix
    indep=-1
    label(n) = 0
    if(els(n)<el(1).or.els(n)>=el(2)) cycle
    do m = 1, n-1
      if(els(n)==els(m)) then
        indep = m
        go to 20
      end if
    end do
20  if(indep==-1) then 
      nmax = nmax + 1
      label(n) = nmax 
    else 
      label(n) = label(indep)
    end if
  end do

end subroutine power_labeling

subroutine calccl_1d(D,alm1,alm2,Cl)
  implicit none
  !I/O
  double precision, intent(in) :: D(2)
  double precision, intent(out) :: Cl(:)
  double complex, intent(in) :: alm1(:), alm2(:)

  Cl = (alm1*conjg(alm2)+alm2*conjg(alm1))*0.5d0/(D(1)*D(2))  ! divided by delta(l=0)

end subroutine calccl_flat_alm1d


subroutine calccl_2d(D,alm1,alm2,Cl)
  implicit none
  !I/O
  double precision, intent(in) :: D(2)
  double precision, intent(out) :: Cl(:)
  double complex, intent(in) :: alm1(:,:), alm2(:,:)
  !internal
  integer :: i, j, n, npix
  double precision, allocatable :: C(:)

  npix = size(alm1)
  n = 1
  do i = 1, size(alm1,dim=1)
    do j = 1, size(alm1,dim=2)
      Cl(n) = (alm1(i,j)*conjg(alm2(i,j))+alm2(i,j)*conjg(alm1(i,j)))*0.5d0/(D(1)*D(2))  ! divided by delta(l=0)
      n = n + 1
    end do
  end do

end subroutine calccl_flat_alm2d

#endif all

end module pstool

