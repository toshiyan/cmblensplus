!////////////////////////////////////////////////////!
! * dft module
!////////////////////////////////////////////////////!

module fftw
  use constants, only: dlc, iu, pi, twopi
  use grid2d,    only: spin_weight
  implicit none

  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)

  interface dft
    module procedure dft_1darray, dft_2darray
  end interface

  interface dft_pol
    module procedure dft_pol_1darray, dft_pol_2darray
  end interface

  !interface derivemap
  !  module procedure derivemap_all, derivemap_nth_1d, derivemap_nth_2d, derivemap_alm_nth_1d
  !end interface

  private FFTW_ESTIMATE
  private dlc, iu, pi, twopi
  private spin_weight

contains 

!////////////////////////////////////////////////////////////////////////!
!
! * FT convention of e.g., Hu & Okamoto (2002)
!
!   E(l) + iB(l) = - \int dx e^(-ilx) [ Q(x) + iU(x) ] e^(-2i\phi)
!   E(l) - iB(l) = - \int dx e^(-ilx) [ Q(x) - iU(x) ] e^(2i\phi)
! 
! Note: dx -> dx/2\pi in Lewis & Challinor (2006) 
!
! * Inversion
!
!   Q(x) + iU(x) = - \int (dl/2pi^2) e^(ilx) [ E(l) + iB(l) ] e^(2i\phi)
!   Q(x) - iU(x) = - \int (dl/2pi^2) e^(ilx) [ E(l) - iB(l) ] e^(-2i\phi)
!
! * Other quantities
!
!   Q(l) + iU(l) = - [ E(l) + iB(l) ] e^(2i\phi) 
!   Q(l) - iU(l) = - [ E(l) - iB(l) ] e^(-2i\phi) 
!
! where Q(l) and U(l) are the FT of Q(x) and U(x), respectively.
! This is equivalent to
! 
!   Q(l) = - E(l)*cos(2\phi) + B(l)*sin(2\phi)  
!   U(l) = - E(l)*sin(2\phi) - B(l)*cos(2\phi)  
!
! or
!
!   E(l) = - Q(l)*cos(2\phi) - U(l)*sin(2\phi)
!   B(l) = Q(l)*sin(2\phi) - U(l)*cos(2\phi)
!
! FT algorithm below assume 
!   
!   Q(x), U(x) <-> Q(l), U(l) <-> E(l), B(l)
!

!////////////////////////////////////////////////////////////////////////!
!
! DESCRIPTION OF dft ROUTINE
!
! [note]
!   trans = -1 -> alm to map
!   trans = +1 -> map to alm
!
! [from FFTW web (http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html) : Sec.4.8.1]
!   An FFTW_FORWARD transform corresponds to a sign of -1 in the exponent of the DFT. 
!   The standard "in-order" output ordering --- the k-th output corresponds to the frequency k/n (or k/T, where T is your total sampling period).
!
!
! * Definition of FT
!
!   FT:  F(lx,ly) = \int dxdy f(x,y)
!   IFT: f(x,y) = \int (dlx/2pi)(dly/2pi) F(lx,ly)
!
! * dx, dl ...
!
!   dx = D(1)/nn(1), dy = D(2)/nn(2)
!   dlx = 2pi/D(1), dly = 2pi/D(2)
!
! * Delta function 
!
!   Delta(lx,ly) = \int dxdy e^{-ixl} 
!       ->  dxdy*n(1)*n(2) = D(1)*D(2)        ( at lx=ly=0 )
!       ->  0                                 ( otherwise )
!
!   Delta(x,y)   = \int (dlx/2pi)*(dly/2pi) e^{ixl} 
!       ->  (dlx/2pi)(dly/2pi)*n(1)*n(2) 
!                   = n(1)*n(2)/D(1)/D(2)     ( at x=y=0 )
!       ->  0                                 ( otherwise )
!
! where the above Delta functions satisfy the usual condition
!
!   \int (dlx/2pi)(dly/2pi) Delta(lx,ly) = 1
!   \int dxdy Delta(x,y)     = 1
!
! The above should be satisfied because 
!
!   \int (dlx/2pi)(dly/2pi) Delta(lx,ly) 
!    = \int (dlx/2pi)(dly/2pi) \int dxdy e^{-ixl}
!    = \int dxdy Delta(x,y)
!

subroutine dft_2darray(map,nn,D,trans)
!* note (Jan 28, 2015)
! Even after FT and inverse FT, the resultant map has tiny error 
! but it biases the angular Cl significantly on large scale, 
! so these values are imposed to zero (add error control). 
!* update (Apr 06, 2019)
! nn(1) and nn(2) can be odd number
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans
  double precision, intent(in) :: D(2)
  complex(dlc), intent(inout) :: map(:,:)
  !internal
  integer :: i, j, m1, m2, n, plan(8)
  double precision :: mean
  complex(dlc) :: f
  complex(dlc), allocatable :: amap(:,:)

  allocate(amap(0:nn(1)-1,0:nn(2)-1))
  do i = 1, nn(1)
    do j = 1, nn(2)
      amap(i-1,j-1) = (-1)**(i+j)*map(i,j)
    end do
  end do

  call DFFTW_PLAN_dft_2D(plan,nn(1),nn(2),amap,amap,trans,FFTW_ESTIMATE)
  call DFFTW_EXEcutE_dft(plan,amap,amap)
  call DFFTW_DESTROY_PLAN(plan)

  ! phase factor
  f = 1d0
  if (mod(nn(1),2)/=0)  f = f*iu
  if (mod(nn(2),2)/=0)  f = f*iu

  ! pre factor
  if (trans==-1) f = conjg(f)/(D(1)*D(2))           ! multiply 1/dS * 1/(nn(1)*nn(2))
  if (trans==1 ) f = f*D(1)*D(2)/dble(nn(1)*nn(2))  ! area of one real-space pixel (dxdy)
 
  m1 = int(nn(1)/2)
  m2 = int(nn(2)/2)
  do i = 1, nn(1)
    do j = 1, nn(2)
      map(i,j) = f*(-1)**(i+j+m1+m2)*amap(i-1,j-1)
    end do
  end do
  deallocate(amap)

  ! error control
  mean = sum(abs(map))/dble(nn(1)*nn(2))
  do i = 1, nn(1)
    do j = 1, nn(2)
      if(abs(map(i,j))<mean*1d-15) map(i,j) = 0d0
    end do
  end do


end subroutine dft_2darray


subroutine dft_1darray(map,nn,D,trans)
!* revised (Jun 8, 2015)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), trans
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(inout) :: map(:)
  !internal
  integer :: i, j, m1, m2, n, plan(8)
  double precision :: mean
  complex(dlc) :: f
  complex(dlc), allocatable :: amap(:,:)

  allocate(amap(0:nn(1)-1,0:nn(2)-1))
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      amap(i-1,j-1) = (-1)**(i+j)*map(n)
      n = n + 1
    end do
  end do

  call DFFTW_PLAN_dft_2D(plan,nn(1),nn(2),amap,amap,trans,FFTW_ESTIMATE)
  call DFFTW_EXEcutE_dft(plan,amap,amap)
  call DFFTW_DESTROY_PLAN(plan)

  ! phase factor
  f = 1d0
  if (mod(nn(1),2)/=0)  f = f*iu
  if (mod(nn(2),2)/=0)  f = f*iu

  ! pre factor
  if (trans==-1) f = conjg(f)/(D(1)*D(2))           ! multiply 1/dS * 1/(nn(1)*nn(2))
  if (trans==1 ) f = f*D(1)*D(2)/dble(nn(1)*nn(2))  ! area of one real-space pixel (dxdy)
 
  m1 = int(nn(1)/2)
  m2 = int(nn(2)/2)
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      map(n) = f*amap(i-1,j-1)*(-1)**(i+j+m1+m2)
      n = n + 1
    end do
  end do
  deallocate(amap)

  ! error control (values smaller than double precision is assigned to zero)
  mean = sum(abs(map))/dble(nn(1)*nn(2))
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      if(abs(map(n))<mean*1d-15) map(n) = 0d0
      n = n + 1
    end do
  end do

end subroutine dft_1darray


subroutine dft_pol_1darray(QU,nn,D,EB,trans)
  !trans=-1 : E,B -> Q+iU
  !trans=+1 : Q+iU -> E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: QU(:,:)
  complex(dlc), intent(inout) :: EB(:,:)
  !internal
  complex(dlc), allocatable :: cP(:), P(:)

  select case(trans)
  case(-1)
    allocate(P(nn(1)*nn(2)))
    P = EB(1,:) + iu*EB(2,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    QU(:,1) = dble(P)
    QU(:,2) = aimag(P)
    deallocate(P)
  case(1)
    allocate(P(nn(1)*nn(2)),cP(nn(1)*nn(2)))
    P = QU(:,1) + iu*QU(:,2)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call dft(cP,nn,D,1)
    call spin_weight(cP,nn,D,1,-2)
    EB(1,:) = (P+cP)*0.5d0
    EB(2,:) = -iu*(P-cP)*0.5d0
    deallocate(P,cP)
  end select

end subroutine dft_pol_1darray


subroutine dft_pol_2darray(QU,nn,D,EB,trans)
  !trans=-1 : E,B -> Q+iU
  !trans=+1 : Q+iU -> E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: QU(:,:,:)
  complex(dlc), intent(inout) :: EB(:,:,:)
  !internal
  complex(dlc), allocatable :: cP(:,:), P(:,:)

  select case(trans)
  case(-1)
    allocate(P(nn(1),nn(2)))
    P = EB(1,:,:) + iu*EB(2,:,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    QU(:,:,1) = dble(P)
    QU(:,:,2) = aimag(P)
    deallocate(P)
  case(1)
    allocate(P(nn(1),nn(2)),cP(nn(1),nn(2)))
    P = QU(:,:,1) + iu*QU(:,:,2)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call dft(cP,nn,D,1)
    call spin_weight(cP,nn,D,1,-2)
    EB(1,:,:) = (P+cP)*0.5d0
    EB(2,:,:) = -iu*(P-cP)*0.5d0
    deallocate(P,cP)
  end select

end subroutine dft_pol_2darray



subroutine dft1d(map0,nx,ny,npix,D,trans,map1)
  implicit none
  integer, intent(in) :: nx, ny, npix, trans
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(npix) :: map0
  double complex, intent(out) :: map1(npix)
  integer :: nn(2)
  double complex :: map(npix)

  nn(1) = nx
  nn(2) = ny
  map   = map0
  call dft_1darray(map,nn,D,trans)
  map1  = map

end subroutine dft1d


subroutine dft2d(map0,nx,ny,D,trans,map1)
  implicit none
  integer, intent(in) :: nx, ny, trans
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: map0
  double complex, intent(out), dimension(nx,ny) :: map1
  integer :: nn(2)
  double complex :: map(nx,ny)

  nn(1) = nx
  nn(2) = ny
  map   = map0
  call dft_2darray(map,nn,D,trans)
  map1  = map

end subroutine dft2d


subroutine dft2dr(map0,nx,ny,D,trans,map1)
  implicit none
  integer, intent(in) :: nx, ny, trans
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: map0
  double precision, intent(out), dimension(nx,ny) :: map1
  integer :: nn(2)
  double complex :: map(nx,ny)

  nn(1) = nx
  nn(2) = ny
  map   = map0
  call dft_2darray(map,nn,D,trans)
  map1  = map

end subroutine dft2dr


subroutine dft2drc(map0,nx,ny,D,trans,map1)
  implicit none
  integer, intent(in) :: nx, ny, trans
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: map0
  double complex, intent(out), dimension(nx,ny) :: map1
  integer :: nn(2)
  double complex :: map(nx,ny)

  nn(1) = nx
  nn(2) = ny
  map   = map0
  call dft_2darray(map,nn,D,trans)
  map1  = map

end subroutine dft2drc


subroutine dft2dcr(map0,nx,ny,D,trans,map1)
  implicit none
  integer, intent(in) :: nx, ny, trans
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: map0
  double precision, intent(out), dimension(nx,ny) :: map1
  integer :: nn(2)
  double complex :: map(nx,ny)

  nn(1) = nx
  nn(2) = ny
  map   = map0
  call dft_2darray(map,nn,D,trans)
  map1  = map

end subroutine dft2dcr



end module fftw




