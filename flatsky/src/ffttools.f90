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
!*  DFT for 1D array. 
!*
!*  Args:
!*    :nx, ny (int)      : Number of x/Lx and y/Ly grids
!*    :npix (int)        : Total number of grids (npix=nx*ny)
!*    :trans (int)       : 1 (map to Fourier) or -1 (Fourier to map)
!*    :D[2] (double)     : Side length (x and y) of map
!*    :map0[pix] (dcmplx): Data on 2D grid to be transformed, with bounds (npix)
!*
!*  Returns:
!*    :map1[pix] (dcmplx): Transformed data on 2D grid, with bounds (npix)
!*
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
!*  DFT for 2D array. 
!*
!*  Args:
!*    :nx, ny (int)      : Number of x/Lx and y/Ly grids
!*    :trans (int)       : 1 (map to Fourier) or -1 (Fourier to map)
!*    :D[2] (double)     : Side length (x and y) of map
!*    :map0[x,y] (dcmplx): Data on 2D grid to be transformed, with bounds (nx,ny)
!*
!*  Returns:
!*    :map1[x,y] (dcmplx): Transformed data on 2D grid, with bounds (nx,ny)
!*
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
!*  DFT for 2D array. 
!*
!*  Args:
!*    :nx, ny (int)      : Number of x/Lx and y/Ly grids
!*    :trans (int)       : 1 (map to Fourier) or -1 (Fourier to map)
!*    :D[2] (double)     : Side length (x and y) of map
!*    :map0[x,y] (double): Data on 2D grid to be transformed, with bounds (nx,ny)
!*
!*  Returns:
!*    :map1[x,y] (double): Transformed data on 2D grid, with bounds (nx,ny)
!*
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
!*  DFT for 2D array. 
!*
!*  Args:
!*    :nx, ny (int)      : Number of x/Lx and y/Ly grids
!*    :trans (int)       : 1 (map to Fourier) or -1 (Fourier to map)
!*    :D[2] (double)     : Side length (x and y) of map
!*    :map0[x,y] (double): Data on 2D grid to be transformed, with bounds (nx,ny)
!*
!*  Returns:
!*    :map1[x,y] (dcmplx): Transformed data on 2D grid, with bounds (nx,ny)
!*
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
!*  DFT for 2D array. 
!*
!*  Args:
!*    :nx, ny (int)      : Number of x/Lx and y/Ly grids
!*    :trans (int)       : 1 (map to Fourier) or -1 (Fourier to map)
!*    :D[2] (double)     : Side length (x and y) of map
!*    :map0[x,y] (dcmplx): Data on 2D grid to be transformed, with bounds (nx,ny)
!*
!*  Returns:
!*    :map1[x,y] (double): Transformed data on 2D grid, with bounds (nx,ny)
!*
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


subroutine dft2dpol(nx,ny,D,Q,U,E,B)
!*  Spin-2 DFT for 2D array. 
!*
!*  Args:
!*    :nx, ny (int)   : Number of x/Lx and y/Ly grids
!*    :D[2] (double)  : Side length (x and y) of map
!*    :Q[x,y] (doub;e): Q on 2D grid to be transformed, with bounds (nx,ny)
!*    :U[x,y] (doub;e): U on 2D grid to be transformed, with bounds (nx,ny)
!*
!*  Returns:
!*    :E[x,y] (dcmplx): E on 2D Fourier grid, with bounds (nx,ny)
!*    :B[x,y] (dcmplx): B on 2D Fourier grid, with bounds (nx,ny)
!*
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: Q, U
  double precision, intent(out), dimension(nx,ny) :: E, B
  integer :: nn(2)
  double precision :: QU(nx,ny,2)
  double complex :: EB(2,nx,ny)

  nn(1) = nx
  nn(2) = ny
  QU(:,:,1) = Q
  QU(:,:,2) = U
  call dft_pol(QU,nn,D,EB,1)
  E  = EB(1,:,:)
  B  = EB(2,:,:)

end subroutine dft2dpol


subroutine idft2dpol(nx,ny,D,E,B,Q,U)
!*  Spin-2 Inverse DFT for 2D array. 
!*
!*  Args:
!*    :nx, ny (int)   : Number of x/Lx and y/Ly grids
!*    :D[2] (double)  : Side length (x and y) of map
!*    :E[x,y] (dcmplx): E on 2D Fourier grid, with bounds (nx,ny)
!*    :B[x,y] (dcmplx): B on 2D Fourier grid, with bounds (nx,ny)
!*
!*  Returns:
!*    :Q[x,y] (doub;e): Q on 2D grid, with bounds (nx,ny)
!*    :U[x,y] (doub;e): U on 2D grid, with bounds (nx,ny)
!*
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: E, B
  double complex, intent(out), dimension(nx,ny) :: Q, U
  integer :: nn(2)
  double precision :: QU(nx,ny,2)
  double complex :: EB(2,nx,ny)

  nn(1) = nx
  nn(2) = ny
  EB(1,:,:) = E
  EB(2,:,:) = B
  call dft_pol(QU,nn,D,EB,-1)
  Q = QU(:,:,1)
  U = QU(:,:,2)

end subroutine idft2dpol


subroutine eb_separate(nx,ny,D,QU,EB,W,Wd)
!*  Compute Smith's pure EB estimator in flatsky
!*
!*  Args:
!*    :nx, ny (int)      : Number of x/Lx and y/Ly grids
!*    :D[2] (double)     : Side length (x and y) of map
!*    :W[x,y] (double)   : Window function, with bounds (nx,ny)
!*    :QU[x,y,2] (double): unmasked Q and U maps, with bounds (nx,ny,2)
!*
!*  Args(optional):
!*    :Wd[5,x,y] (double): Precomputed window function derivaives, dW/dx, dW/dw, d^2W/dx^2, d^2W/dxdy, d^2W/dy^2, with bounds (5,nx,ny)
!*
!*  Returns:
!*    :EB[2,x,y] (dcmplx): E and B modes in 2D Fourier grid, with bounds (2,nx,ny)
!*
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
  !docstr :: Wd = np.zeros((5,nx,ny))
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
  if (present(Wd).and.sum(Wd)/=0d0) then
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

  !Transform to Fourier Space
  call dft(E1,nn,D,1)
  call dft(B1,nn,D,1)
  call dft(E2x,nn,D,1)
  call dft(E2y,nn,D,1)
  call dft(B2x,nn,D,1)
  call dft(B2y,nn,D,1)

  !Add corrections
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


