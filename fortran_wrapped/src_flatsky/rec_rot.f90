!///////////////////////////////////////////////////////////!
! Rotation Reconstruction Kernel
!///////////////////////////////////////////////////////////!

module rec_rot
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qte(nx,ny,D,rL,fC,T,E,alm)
!*  Reconstructing anisotropic pol. rot. angles from TE quadratic estimator
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :fC[lx,ly] (double) : TE cross power spectrum on 2D grid, with bounds (nx,ny)
!*    :T[lx,ly] (dcmplx)  : 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
!*    :E[lx,ly] (dcmplx)  : 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
!*
!*  Returns:
!*    :alm[lx,ly] (dcmplx) : 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, fC, T, E
  !f2py intent(out) alm
  !f2py depend(nx) fC, T, E, alm
  !f2py depend(ny) fC, T, E, alm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T, E
  double complex, intent(out), dimension(nx,ny) :: alm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: aT, aE, ei2p, xlm

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  aT = lmask*fC*T*ei2p
  aE = lmask*E*conjg(ei2p)

  ! convolution
  call dft(aE,nn,D,-1)
  call dft(aT,nn,D,-1)
  xlm = -2d0*aimag(aE*aT)
  call dft(xlm,nn,D,1)

  alm = xlm

end subroutine qte

subroutine qtb(nx,ny,D,rL,fC,T,B,alm)
!*  Reconstructing anisotropic pol. rot. angles from TB quadratic estimator
!*
!*  Args:
!*    :nx, ny (int)         : Number of Lx and Ly grids
!*    :D[xy] (double)       : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)          : Minimum and maximum multipole of CMB for reconstruction
!*    :fC[lx,ly] (double)   : TE cross power spectrum on 2D grid, with bounds (nx,ny)
!*    :T[lx,ly] (dcmplx)    : 2D Fourier modes of inverse-variance filtered temperature, with bounds (nx,ny)
!*    :B[lx,ly] (dcmplx)    : 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)
!*
!*  Returns:
!*    :alm[lx,ly] (dcmplx) : 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, fC, T, B
  !f2py intent(out) alm
  !f2py depend(nx) fC, T, B, alm
  !f2py depend(ny) fC, T, B, alm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T, B
  double complex, intent(out), dimension(nx,ny) :: alm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wB, wT, xlm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wT = lmask*fC*T*ei2p
  wB = lmask*B*conjg(ei2p)

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wT,nn,D,-1)
  xlm = 2d0*dble(wB*wT)
  call dft(xlm,nn,D,1)

  alm = xlm

end subroutine qtb

subroutine qee(nx,ny,D,rL,fC,E1,E2,alm)
!*  Reconstructing anisotropic pol. rot. angles from EE quadratic estimator
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :fC[lx,ly] (double) : EE power spectrum on 2D grid, with bounds (nx,ny)
!*    :E1[lx,ly] (dcmplx) : 2D Fourier modes of 1st inverse-variance filtered E-mode, with bounds (nx,ny)
!*    :E2[lx,ly] (dcmplx) : 2D Fourier modes of 2nd inverse-variance filtered E-mode, with bounds (nx,ny)
!*
!*  Returns:
!*    :alm[lx,ly] (dcmplx): 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, fC, E1, E2
  !f2py intent(out) alm
  !f2py depend(nx) fC, E1, E2, alm
  !f2py depend(ny) fC, E1, E2, alm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: E1, E2
  double complex, intent(out), dimension(nx,ny) :: alm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wE1, wE2, xlm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wE1 = lmask*E1*conjg(ei2p)
  wE2 = lmask*fC*E2*ei2p

  ! convolution
  call dft(wE1,nn,D,-1)
  call dft(wE2,nn,D,-1)
  xlm = -2d0*aimag(wE1*wE2)
  call dft(xlm,nn,D,1)

  alm = xlm

end subroutine qee

subroutine qeb(nx,ny,D,rL,EE,E,B,alm,BB)
!*  Reconstructing anisotropic pol. rot. angles from EB quadratic estimator
!*
!*  Args:
!*    :nx, ny (int)         : Number of Lx and Ly grids
!*    :D[xy] (double)       : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)          : Minimum and maximum multipole of CMB for reconstruction
!*    :fC[lx,ly] (double)   : EE power spectrum on 2D grid, with bounds (nx,ny)
!*    :E[lx,ly] (dcmplx)    : 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
!*    :B[lx,ly] (dcmplx)    : 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)
!*
!*  Args(Optional):
!*    :BB[lx,ly] (double) : Theory B-mode spectrum on 2D grid, with bounds (nx,ny), default to BB=0
!*
!*  Returns:
!*    :alm[lx,ly] (dcmplx) : 2D Fourier modes of anisotropic pol. rot. angles, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, EE, E, B, BB
  !f2py intent(out) alm
  !f2py depend(nx) EE, E, B, alm, BB
  !f2py depend(ny) EE, E, B, alm, BB
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: EE
  double complex, intent(in), dimension(nx,ny) :: E, B
  double complex, intent(out), dimension(nx,ny) :: alm
  !(optional)
  double precision, intent(in), optional, dimension(nx,ny) :: BB
  !opt4py :: BB = 0
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: lmask, els
  double complex, dimension(nx,ny) :: wE, wB, alm1, alm2, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  !filtering
  call make_lmask(nn,D,rL,lmask)
  wE = lmask*EE*E*ei2p
  wB = lmask*B*conjg(ei2p)

  !convolution
  call dft(wE,nn,D,-1)
  call dft(wB,nn,D,-1)
  alm1 = 2d0*dble(wE*wB)
  call dft(alm1,nn,D,1)

  if (sum(abs(BB))/=0d0) then
    !filtering
    wE = lmask*E*ei2p
    wB = lmask*BB*B*conjg(ei2p)
    !convolution
    call dft(wE,nn,D,-1)
    call dft(wB,nn,D,-1)
    alm2 = 2d0*dble(wE*wB)
    call dft(alm2,nn,D,1)
  end if

  !estimator
  alm = alm1 + alm2

end subroutine qeb



end module rec_rot



