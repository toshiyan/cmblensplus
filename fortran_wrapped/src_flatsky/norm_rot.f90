!///////////////////////////////////////////////////////////!
! Normalization of pol. rot. reconstruction
!///////////////////////////////////////////////////////////!

module norm_rot
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qeb(nx,ny,D,rL,IE,IB,EE,eL,Aa,BB)
!*  Normalization of the EB quadratic estimator for anisotropic pol. rot. angles
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
!*  Args(Optional):
!*    :BB[lx,ly] (double) : Theory B-mode spectrum on 2D grid, with bounds (nx,ny), default to BB=0
!*
!*  Returns:
!*    :Aa[lx,ly] (dcmplx) : Normalization of anisotropic pol. rot. angles on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, rL, eL, D, EE, IE, IB, BB
  !f2py intent(out) Aa
  !f2py depend(nx) EE, IE, IB, Aa, BB
  !f2py depend(ny) EE, IE, IB, Aa, BB
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL, eL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: EE, IE, IB
  double precision, intent(out), dimension(nx,ny) :: Aa
  !(optional)
  double precision, intent(in), optional, dimension(nx,ny) :: BB
  !f2py double precision :: BB = 0
  !internal
  integer :: i, j, nn(2)
  double precision, dimension(nx,ny) :: els, lmask, Aaa
  double complex, dimension(nx,ny) :: ei2p
  double complex, dimension(nx,ny) :: Al, Bl, Cl, X1, X2, Y1, Y2

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  call make_lmask(nn,D,rL,lmask)

  !EE^2 part
  X1 = lmask*IE*EE**2
  X2 = lmask*IE*EE**2*ei2p**2
  Y1 = lmask*IB
  Y2 = lmask*IB*conjg(ei2p)**2
  call dft(X1,nn,D,-1)
  call dft(X2,nn,D,-1)
  call dft(Y1,nn,D,-1)
  call dft(Y2,nn,D,-1)
  Al = X1*Y1 + X2*Y2
  call dft(Al,nn,D,1)
  Al = dble(Al)

  if (present(BB).and.sum(abs(BB))/=0d0) then
    !BB^2 part
    X1 = lmask*IE
    X2 = lmask*IE*ei2p**2
    Y1 = lmask*IB*BB**2
    Y2 = lmask*IB*BB**2*conjg(ei2p)**2
    call dft(X1,nn,D,-1)
    call dft(X2,nn,D,-1)
    call dft(Y1,nn,D,-1)
    call dft(Y2,nn,D,-1)
    Bl = X1*Y1 + X2*Y2
    call dft(Bl,nn,D,1)
    Bl = dble(Bl)
    !EExBB part
    X1 = lmask*IE*EE
    X2 = lmask*IE*EE*ei2p**2
    Y1 = lmask*IB*BB
    Y2 = lmask*IB*BB*conjg(ei2p)**2
    call dft(X1,nn,D,-1)
    call dft(X2,nn,D,-1)
    call dft(Y1,nn,D,-1)
    call dft(Y2,nn,D,-1)
    Cl = X1*Y1 + X2*Y2
    call dft(Cl,nn,D,1)
    Cl = dble(Cl)
  end if

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Aaa = 2d0*(Al+Bl) - 4d0*Cl

  ! inversion
  Aa = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Aaa(i,j)>0)  Aa(i,j) = 1d0/Aaa(i,j)
    end do
  end do


end subroutine qeb



end module norm_rot



