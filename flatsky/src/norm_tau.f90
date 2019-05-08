!///////////////////////////////////////////////////////////!
! Normalization of tau teconstruction
!///////////////////////////////////////////////////////////!

module norm_tau
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qeb(nx,ny,D,rL,IE,IB,EE,eL,At)
!*  Normalization of the EB quadratic estimator for patchy tau
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :IE[lx,ly] (double) : Inverse of observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :IB[lx,ly] (double) : Inverse of observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :EE[lx,ly] (double) : Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :At[lx,ly] (dcmplx) : Normalization of patchy tau on 2D grid, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: IE, IB, EE
  double precision, intent(out), dimension(nx,ny) :: At
  !internal
  integer :: i, j, nn(2)
  double precision, dimension(nx,ny) :: lmask, Att, els
  double complex, dimension(nx,ny) :: X1, X2, Y1, Y2, Al, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  !filtering
  X1 = lmask*EE**2*IE
  X2 = lmask*EE**2*IE*ei2p**2
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
  Att = Al*0.5d0

  !inversion
  At = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Att(i,j)>0)  At(i,j) = 1d0/Att(i,j)
    end do
  end do

end subroutine qeb


end module norm_tau


