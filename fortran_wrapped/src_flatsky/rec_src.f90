!///////////////////////////////////////////////////////////!
! Point source reconstruction kernel
!///////////////////////////////////////////////////////////!

module rec_src
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qtt(nx,ny,D,rL,T1,T2,slm)
!*  Reconstructing point source fields from the temperature quadratic estimator
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :fC[lx,ly] (double) : Temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :T1[lx,ly] (dcmplx) : 2D Fourier modes of 1st inverse-variance filtered temperature, with bounds (nx,ny)
!*    :T2[lx,ly] (dcmplx) : 2D Fourier modes of 2nd inverse-variance filtered temperature, with bounds (nx,ny)
!*
!*  Returns:
!*    :slm[lx,ly] (dcmplx): 2D Fourier modes of point source fields, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, T1, T2
  !f2py intent(out) slm
  !f2py depend(nx) T1, T2, slm
  !f2py depend(ny) T1, T2, slm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: T1, T2
  double complex, intent(out), dimension(nx,ny) :: slm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: aT, alm

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els)

  !filtering
  call make_lmask(nn,D,rL,lmask)
  aT  = lmask*T1
  alm = lmask*T2

  !convolution
  call dft(aT,nn,D,-1)
  call dft(alm,nn,D,-1)
  alm = aT*alm
  call dft(alm,nn,D,1)

  slm = alm

end subroutine qtt

subroutine qte(nx,ny,D,rL,T,E,slm)
!*  Reconstructing point source fields from the suboptimal TE quadratic estimator
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
!*    :slm[lx,ly] (dcmplx) : 2D Fourier modes of point source fields, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, T, E
  !f2py intent(out) slm
  !f2py depend(nx) T, E, slm
  !f2py depend(ny) T, E, slm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: T, E
  double complex, intent(out), dimension(nx,ny) :: slm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: aT, aE, bT, bE, ei2p, alm

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  aE = lmask*E*conjg(ei2p)
  bE = lmask*T*ei2p
  aT = lmask*T
  bT = lmask*E

  ! convolution
  call dft(aE,nn,D,-1)
  call dft(aT,nn,D,-1)
  call dft(bE,nn,D,-1)
  call dft(bT,nn,D,-1)
  alm = aimag(aE*bE) - aT*bT
  call dft(alm,nn,D,1)

  slm = alm

end subroutine qte

subroutine qtb(nx,ny,D,rL,T,B,slm)
!*  Reconstructing point source fields from the TB quadratic estimator
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
!*    :slm[lx,ly] (dcmplx) : 2D Fourier modes of point source fields, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, T, B
  !f2py intent(out) slm
  !f2py depend(nx) T, B, slm
  !f2py depend(ny) T, B, slm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: T, B
  double complex, intent(out), dimension(nx,ny) :: slm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wB, wT, alm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wB = lmask*B*conjg(ei2p)
  wT = lmask*T*ei2p

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wT,nn,D,-1)
  alm = aimag(wB*wT)
  call dft(alm,nn,D,1)

  slm = alm

end subroutine qtb

subroutine qee(nx,ny,D,rL,E1,E2,slm)
!*  Reconstructing point source fields from the EE quadratic estimator
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
!*    :slm[lx,ly] (dcmplx) : 2D Fourier modes of point source fields, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, E1, E2
  !f2py intent(out) slm
  !f2py depend(nx) E1, E2, slm
  !f2py depend(ny) E1, E2, slm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: E1, E2
  double complex, intent(out), dimension(nx,ny) :: slm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wE1, wE2, alm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wE1 = lmask*E1*conjg(ei2p)
  wE2 = lmask*E2*ei2p

  ! convolution
  call dft(wE1,nn,D,-1)
  call dft(wE2,nn,D,-1)
  alm = real(wE1*wE2)
  call dft(alm,nn,D,1)

  slm = alm

end subroutine qee

subroutine qeb(nx,ny,D,rL,E,B,slm)
!*  Reconstructing point source fields from the EB quadratic estimator
!*
!*  Args:
!*    :nx, ny (int)         : Number of Lx and Ly grids
!*    :D[xy] (double)       : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)          : Minimum and maximum multipole of CMB for reconstruction
!*    :fC[lx,ly] (double)   : EE power spectrum on 2D grid, with bounds (nx,ny)
!*    :E[lx,ly] (dcmplx)    : 2D Fourier modes of inverse-variance filtered E-mode, with bounds (nx,ny)
!*    :B[lx,ly] (dcmplx)    : 2D Fourier modes of inverse-variance filtered B-mode, with bounds (nx,ny)
!*
!*  Returns:
!*    :slm[lx,ly] (dcmplx) : 2D Fourier modes of point source fields, with bounds (nx,ny)
!*
  implicit none
  !f2py intent(in) nx, ny, rL, D, E, B
  !f2py intent(out) slm
  !f2py depend(nx) E, B, slm
  !f2py depend(ny) E, B, slm
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: E, B
  double complex, intent(out), dimension(nx,ny) :: slm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wE, wB, alm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wE = lmask*E*ei2p
  wB = lmask*B*conjg(ei2p)

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wE,nn,D,-1)
  alm = aimag(wB*wE)
  call dft(alm,nn,D,1)

  slm = alm

end subroutine qeb



end module rec_src



