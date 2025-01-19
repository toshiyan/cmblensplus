!///////////////////////////////////////////////////////////!
! Normalization of tau teconstruction
!///////////////////////////////////////////////////////////!

module norm_src
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qtt(nx,ny,D,rL,OT,eL,As)
!*  Normalization of the temperature quadratic estimator for point source
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OT[lx,ly] (double) : Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :As[lx,ly] (dcmplx) : Normalization of point source on 2D grid, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: OT
  double precision, intent(out), dimension(nx,ny) :: As
  !internal
  integer :: i, j, nn(2)
  double precision, dimension(nx,ny) :: lmask, lx, ly, Al
  double complex, dimension(nx,ny) :: A

  nn = (/nx,ny/)

  call elarrays_2d(nn,D,elx=lx,ely=ly)

  ! filtering
  call make_lmask(nn,D,rL,lmask)

  A  = lmask * OT

  ! convolution
  call dft(A,nn,D,-1)
  A = A**2
  call dft(A,nn,D,1)

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Al = lmask * A

  ! inversion
  As = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Al(i,j)>0)  As(i,j) = 1d0/Al(i,j)
    end do
  end do

end subroutine qtt


subroutine qeb(nx,ny,D,rL,IE,IB,eL,As)
!*  Normalization of the EB quadratic estimator for point source 
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :IE[lx,ly] (double) : Inverse of observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :IB[lx,ly] (double) : Inverse of observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :As[lx,ly] (dcmplx) : Normalization of point source on 2D grid, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: IE, IB
  double precision, intent(out), dimension(nx,ny) :: As
  !internal
  integer :: i, j, nn(2)
  double precision, dimension(nx,ny) :: lmask, Ass, els
  double complex, dimension(nx,ny) :: X1, X2, Y1, Y2, Al, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  !filtering
  X1 = lmask*IE
  X2 = lmask*IE*ei2p**2
  Y1 = lmask*IB
  Y2 = lmask*IB*conjg(ei2p**2)

  !convolution
  call dft(X1,nn,D,-1)
  call dft(X2,nn,D,-1)
  call dft(Y1,nn,D,-1)
  call dft(Y2,nn,D,-1)
  Al = X1*Y1 - dble(X2*Y2)
  call dft(Al,nn,D,1)

  !normalization
  call make_lmask(nn,D,eL,lmask)
  Ass = lmask * Al*0.5d0

  !inversion
  As = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Ass(i,j)>0)  As(i,j) = 1d0/Ass(i,j)
    end do
  end do

end subroutine qeb


end module norm_src


