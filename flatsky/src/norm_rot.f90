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
!*  Return normalization of EB quadratic estimator for anisotropic pol. rot. angles
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL, eL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: EE, IE, IB
  double precision, intent(out), dimension(nx,ny) :: Aa
  !(optional)
  double precision, intent(in), optional, dimension(nx,ny) :: BB
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


