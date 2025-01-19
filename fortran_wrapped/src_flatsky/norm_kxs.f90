!///////////////////////////////////////////////////////////!
! 2D Cross Normalization
!///////////////////////////////////////////////////////////!

module norm_kxs
  use constants, only: iu, i1, pi, dlc, twopi
  use grid2d, only: elarrays_2d, make_lmask
  use fftw,  only: dft
  implicit none

  private iu, i1, pi, dlc, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qtt(nx,ny,D,rL,OT,TT,eL,Aks)
!*  Normalization of the temperature quadratic estimators between CMB lensing potential and patchy tau
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OT[lx,ly] (double) : Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :TT[lx,ly] (double) : Theory temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :Aks[lx,ly] (dcmplx) : Lensing-tau cross normalization on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, eL, rL, D, OT, TT
  !f2py intent(out) Aks
  !f2py depend(nx) OT, TT, Aks
  !f2py depend(ny) OT, TT, Aks
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: OT, TT
  double precision, intent(out), dimension(nx,ny) :: Aks
  !internal
  integer :: i, j, nn(2)
  double precision, dimension(nx,ny) :: lmask, lx, ly, Alks
  double complex, dimension(nx,ny) :: A, Ax, Ay

  nn = (/nx,ny/)

  call elarrays_2d(nn,D,elx=lx,ely=ly)

  ! filtering
  call make_lmask(nn,D,rL,lmask)

  Ax = lmask * lx * TT*OT
  Ay = lmask * ly * TT*OT
  A  = lmask * OT

  ! convolution
  call dft(Ax,nn,D,-1)
  call dft(Ay,nn,D,-1)
  call dft(A,nn,D,-1)
  Ax = A*Ax
  Ay = A*Ay

  call dft(Ax,nn,D,1)
  call dft(Ay,nn,D,1)

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Alks = lmask * (lx*Ax + ly*Ay)

  ! inversion
  Aks = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Alks(i,j)>0)  Aks(i,j) = 1d0/Alks(i,j)
    end do
  end do

end subroutine qtt



end module norm_kxs



