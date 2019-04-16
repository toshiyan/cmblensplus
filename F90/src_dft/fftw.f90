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

  interface derivemap
    module procedure derivemap_all, derivemap_nth_1d, derivemap_nth_2d, derivemap_alm_nth_1d
  end interface

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
  complex(dlc) :: amap(0:nn(1)-1,0:nn(2)-1), f


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


subroutine window_deriv(nx,ny,D,W,Wx,Wy,Wxx,Wxy,Wyy)
! for apodization window
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: W
  double complex, intent(inout), dimension(nx,ny) :: Wx, Wy, Wxx, Wxy, Wyy
  integer :: i, j, n, npix, nn(2)
  double precision :: lx, ly
  double complex, allocatable :: Wl(:,:)

  nn   =(/nx,ny/)
  npix = nn(1)*nn(2)

  !* array in Fourier space
  allocate(Wl(nx,ny))
  Wl = W
  call dft(Wl,nn,D,1)
  do i = 1, nx
    lx = twopi*dble(i-1-nx*0.5d0)/D(1)
    do j = 1, ny
      ly = twopi*dble(j-1-ny*0.5d0)/D(2)
      Wx(i,j)  = iu*lx*Wl(i,j)
      Wy(i,j)  = iu*ly*Wl(i,j)
      Wxx(i,j) = -lx**2*Wl(i,j)
      Wxy(i,j) = -lx*ly*Wl(i,j)
      Wyy(i,j) = -ly**2*Wl(i,j)
    end do
  end do 
  deallocate(Wl)

  !* Derivatives of windows functions
  call dft(Wx,nn,D,-1)
  call dft(Wy,nn,D,-1)
  call dft(Wxx,nn,D,-1)
  call dft(Wxy,nn,D,-1)
  call dft(Wyy,nn,D,-1)

end subroutine window_deriv


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


subroutine eb_separate(nx,ny,D,QU,EB,W,Wd)
! * compute Smith's pure EB estimator
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: W
  double precision, intent(in), dimension(nx,ny,2) :: QU
  double complex, intent(out), dimension(2,nx,ny) :: EB
  double precision, intent(in), optional, dimension(5,nx,ny) :: Wd
  !internal
  integer :: i, j, n, npix, nn(2)
  double precision :: lx, ly
  double precision, allocatable :: WP(:,:,:)
  double complex, allocatable, dimension(:,:) :: Wl,E1,E2x,E2y,B1,B2x,B2y
  double complex, allocatable, dimension(:,:) :: Wx,Wy,Wxx,Wxy,Wyy
  double complex, allocatable, dimension(:,:,:) :: WEB

  nn   = (/nx,ny/)
  npix = nx*ny

  allocate(WP(nx,ny,2),WEB(2,nx,ny))
  WP(:,:,1) = QU(:,:,1)*W
  WP(:,:,2) = QU(:,:,2)*W
  call dft_pol(WP,nn,D,WEB,1)
  deallocate(WP)

  !//// derivatives of window functions ////!
  allocate(Wx(nx,ny),Wy(nx,ny),Wxx(nx,ny),Wxy(nx,ny),Wyy(nx,ny))
  if (present(Wd)) then
    Wx  = Wd(1,:,:)
    Wy  = Wd(2,:,:)
    Wxx = Wd(3,:,:)
    Wyy = Wd(4,:,:)
    Wxy = Wd(5,:,:)
  else
    call window_deriv(nx,ny,D,W,Wx,Wy,Wxx,Wxy,Wyy)
  end if 

  !//// correction terms ////!
  allocate(E1(nx,ny),B1(nx,ny),E2x(nx,ny),E2y(nx,ny),B2x(nx,ny),B2y(nx,ny))
  E1  = QU(:,:,1)*(Wyy-Wxx) - 2*QU(:,:,2)*Wxy
  B1  = QU(:,:,2)*(Wyy-Wxx) + 2*QU(:,:,1)*Wxy
  E2x =  2*iu*( QU(:,:,1)*Wx+QU(:,:,2)*Wy )
  E2y =  2*iu*( QU(:,:,2)*Wx-QU(:,:,1)*Wy )
  B2x =  2*iu*( QU(:,:,2)*Wx-QU(:,:,1)*Wy ) 
  B2y = -2*iu*( QU(:,:,1)*Wx+QU(:,:,2)*Wy ) 
  deallocate(Wx,Wy,Wxx,Wxy,Wyy)

  !* Transform to Fourier Space
  call dft(E1,nn,D,1)
  call dft(B1,nn,D,1)
  call dft(E2x,nn,D,1)
  call dft(E2y,nn,D,1)
  call dft(B2x,nn,D,1)
  call dft(B2y,nn,D,1)

  !* add corrections
  EB = 0d0
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(lx==0.and.ly==0) cycle
      EB(1,i,j) = WEB(1,i,j) + (E1(i,j)+lx*E2x(i,j)+ly*E2y(i,j))/(lx**2+ly**2)
      EB(2,i,j) = WEB(2,i,j) + (B1(i,j)+lx*B2x(i,j)+ly*B2y(i,j))/(lx**2+ly**2)
    end do
  end do
  deallocate(E1,B1,E2x,E2y,B2x,B2y,WEB)

end subroutine eb_separate


!//// nth order derivatives ////!

subroutine derivemap_nth_1d(nn,D,map,dmap,nth)
! return nth order derivatives of input map
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nth
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(:) :: map
  double precision, intent(out) :: dmap(:,:)
  !internal
  integer :: i, j, k, n, npix
  double precision :: lx, ly
  double complex, allocatable :: alm(:), dftmap(:)

  npix = nn(1)*nn(2)

  allocate(alm(npix))
  alm = cmplx(map)
  call dft(alm,nn,D,1)

  ! nth order (n>=1)
  do k = 0, nth
    n = 1
    allocate(dftmap(npix))
    do i = 1, nn(1)
      lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      do j = 1, nn(2)
        ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
        dftmap(n) = (-iu)**nth*lx**(nth-k)*ly**k*alm(n) !need minus sign of iu for code consistency
        n = n + 1
      end do
    end do
    ! inverse-FT 
    call dft(dftmap,nn,D,-1)
    dmap(k+1,:) = dftmap
    deallocate(dftmap)
  end do

  deallocate(alm)

end subroutine derivemap_nth_1d


subroutine derivemap_nth_2d(D,map,derivmap,nth)
  implicit none
  !I/O
  integer, intent(in) :: nth
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(:,:) :: map
  double precision, intent(out) :: derivmap(:,:,:)
  !internal
  integer :: i,j,k,nn(2)
  double precision :: lx,ly
  double complex, allocatable :: ftmap(:,:), dftmap(:,:)

  ! get information
  nn(1) = size(map,dim=1)
  nn(2) = size(map,dim=2)

  ! prepare FT map
  allocate(ftmap(nn(1),nn(2))); ftmap=0d0
  ftmap = map
  call dft(ftmap,nn,D,1)

  ! nth order (n>=1)
  do k = 0, nth

    ! compute higher-order derivatives in F-space
    allocate(dftmap(nn(1),nn(2)))
    do i = 1, nn(1)
      lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      do j = 1, nn(2)
        ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
        dftmap(i,j) = (-iu)**nth*lx**(nth-k)*ly**k*ftmap(i,j) !need minus sign of iu for code consistency
      end do
    end do

    ! inverse-FT 
    call dft(dftmap,nn,D,-1)
    derivmap(k+1,:,:) = dftmap
    deallocate(dftmap)

  end do

end subroutine derivemap_nth_2d


subroutine derivemap_alm_nth_1d(nn,D,alm,dmap,nth)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nth
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in) :: alm(:)
  double precision, intent(out) :: dmap(:,:)
  !internal
  integer :: i, j, k, n, npix
  double precision :: lx, ly
  double complex, allocatable :: dftmap(:)

  npix = nn(1)*nn(2)

  ! nth order (n>=1)
  do k = 0, nth
    n = 1
    allocate(dftmap(npix))
    do i = 1, nn(1)
      lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      do j = 1, nn(2)
        ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
        dftmap(n) = (-iu)**nth*lx**(nth-k)*ly**k*alm(n) !need minus sign of iu for code consistency
        n = n + 1
      end do
    end do
    ! inverse-FT 
    call dft(dftmap,nn,D,-1)
    dmap(k+1,:) = dftmap
    deallocate(dftmap)
  end do

end subroutine derivemap_alm_nth_1d


subroutine derivemap_all(D,map,derivmap)
! return derivatives up to nth order
  implicit none
  !I/O
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(:,:) :: map
  double precision, intent(out) :: derivmap(:,:,:,:)
  !internal
  integer :: i,j,k,n,nth,nn(2)
  double precision :: lx,ly
  double complex, allocatable :: ftmap(:,:), dftmap(:,:)

  ! get information
  nn(1) = size(map,dim=1)
  nn(2) = size(map,dim=2)
  nth   = size(derivmap,dim=1) - 1

  ! prepare FT map
  allocate(ftmap(nn(1),nn(2))); ftmap=0d0
  ftmap = map
  call dft(ftmap,nn,D,1)

  ! 0th order map
  derivmap(1,1,:,:) = map

  ! nth order (n>=1)
  do n = 1, nth
    do k = 0, n

      ! compute higher-order derivatives in F-space
      allocate(dftmap(nn(1),nn(2)))
      do i = 1, nn(1)
        lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
        do j = 1, nn(2)
          ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
          dftmap(i,j) = (-iu)**n*lx**(n-k)*ly**k*ftmap(i,j) !need minus sign of iu for code consistency
        end do
      end do

      ! inverse-FT 
      call dft(dftmap,nn,D,-1)
      derivmap(k+1,n+1,:,:) = dftmap
      deallocate(dftmap)

    end do
  end do

end subroutine derivemap_all


end module fftw




