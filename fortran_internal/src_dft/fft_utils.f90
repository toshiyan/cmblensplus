!////////////////////////////////////////////////////!
! * dft utilities
!////////////////////////////////////////////////////!

module fft_utils
  use constants, only: dlc, iu, pi, twopi
  use grid2d,    only: make_binmask
  use fftw,      only: dft, dft_pol
  implicit none

  interface derivemap
    module procedure derivemap_all, derivemap_nth_1d, derivemap_nth_2d, derivemap_alm_nth_1d
  end interface

  private dlc, iu, pi, twopi
  private make_binmask
  private dft, dft_pol

contains 


subroutine eb_separate(nx,ny,D,QU,EB,W,Wd)
! * compute the chi estimator
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



subroutine prep_filtered_map(nx,ny,D,bp,alm,fmap,fl0)
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  double precision, intent(in) :: D(2), bp(:)
  double complex, intent(in) :: alm(:,:,:)
  double precision, intent(out) :: fmap(:,:,:,:)
  !optional
  double precision, intent(in), optional :: fl0(:,:,:)
  !internal
  integer :: b, s, bn, smax, nn(2)
  double precision, allocatable :: fl(:,:,:)
  double complex, allocatable :: falm(:,:)

  bn   = size(bp) - 1
  smax = size(alm,dim=1)
  nn   = (/nx,ny/)

  ! create bin filter
  allocate(fl(bn,nx,ny)); fl=0d0
  if (present(fl0)) then
    fl = fl0
  else
    call make_binmask(nn,D,bp(1:bn),bp(2:bn+1),fl,bn)
  end if

  ! compute filtered map
  fmap = 0d0
  do b = 1, bn
    do s = 1, smax
      allocate(falm(nx,ny)); falm=0d0
      falm = fl(b,:,:)*alm(s,:,:)
      call dft(falm,nn,D,-1)
      fmap(s,b,:,:) = dble(falm)
      deallocate(falm)
    end do
  end do

  deallocate(fl)

end subroutine prep_filtered_map





end module fft_utils