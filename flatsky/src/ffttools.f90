!////////////////////////////////////////////////////!
! * dft module
!////////////////////////////////////////////////////!

module ffttools
  use constants, only: iu, pi, twopi
  use fftw,      only: dft, dft_pol, window_deriv
  implicit none

  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)

  interface dft
    module procedure :: dft1d, dft2d, dft2dr, dft2drc, dft2dcr
  end interface dft

  private FFTW_ESTIMATE
  private iu, pi, twopi
  private dft, dft_pol, window_deriv

contains 


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
  call dft(map,nn,D,trans)
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
  call dft(map,nn,D,trans)
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
  call dft(map,nn,D,trans)
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
  call dft(map,nn,D,trans)
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
  call dft(map,nn,D,trans)
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
  !(optional)
  double precision, intent(in), optional, dimension(5,nx,ny) :: Wd
  !f2py double precision :: Wd = 0
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


end module ffttools


