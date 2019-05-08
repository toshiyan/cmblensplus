!///////////////////////////////////////////////////////////!
! Tau reconstruction kernel
!///////////////////////////////////////////////////////////!

module rec_tau
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qtt(nx,ny,D,rL,fC,T1,T2,tlm)
!*  Reconstructing patchy tau from the temperature quadratic estimator
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
!*    :tlm[lx,ly] (dcmplx): 2D Fourier modes of patchy tau, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T1, T2
  double complex, intent(out), dimension(nx,ny) :: tlm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: aT, alm

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els)

  !filtering
  call make_lmask(nn,D,rL,lmask)
  aT  = lmask*T1
  alm = lmask*fC*T2

  !convolution
  call dft(aT,nn,D,-1)
  call dft(alm,nn,D,-1)
  alm = aT*alm
  call dft(alm,nn,D,1)

  tlm = alm

end subroutine qtt


subroutine qte(nx,ny,D,rL,fC,T,E,tlm)
!*  Reconstructing patchy tau from the suboptimal TE quadratic estimator
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
!*    :tlm[lx,ly] (dcmplx) : 2D Fourier modes of patchy tau, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T, E
  double complex, intent(out), dimension(nx,ny) :: tlm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: aT, aE, bT, bE, ei2p, alm

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  aT = lmask*fC*T*ei2p
  aE = lmask*E*conjg(ei2p)
  bT = lmask*T
  bE = lmask*fC*E

  ! convolution
  call dft(aT,nn,D,-1)
  call dft(aE,nn,D,-1)
  call dft(bT,nn,D,-1)
  call dft(bE,nn,D,-1)
  alm = dble(aT*aE) + bT*bE
  call dft(alm,nn,D,1)

  tlm = alm

end subroutine qte


subroutine qtb(nx,ny,D,rL,fC,T,B,tlm)
!*  Reconstructing patchy tau from the TB quadratic estimator
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
!*    :tlm[lx,ly] (dcmplx) : 2D Fourier modes of patchy tau, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T, B
  double complex, intent(out), dimension(nx,ny) :: tlm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wB, wT, alm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wB = lmask*B*conjg(ei2p)
  wT = lmask*fC*T*ei2p

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wT,nn,D,-1)
  alm = aimag(wB*wT)
  call dft(alm,nn,D,1)

  tlm = alm

end subroutine qtb


subroutine qee(nx,ny,D,rL,fC,E1,E2,tlm)
!*  Reconstructing patchy tau from the EE quadratic estimator
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
!*    :tlm[lx,ly] (dcmplx) : 2D Fourier modes of patchy tau, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: E1, E2
  double complex, intent(out), dimension(nx,ny) :: tlm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wE1, wE2, alm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wE1 = lmask*E1*conjg(ei2p)
  wE2 = lmask*fC*E2*ei2p

  ! convolution
  call dft(wE1,nn,D,-1)
  call dft(wE2,nn,D,-1)
  alm = dble(wE1*wE2)
  call dft(alm,nn,D,1)

  tlm = alm

end subroutine qee


subroutine qeb(nx,ny,D,rL,fC,E,B,tlm)
!*  Reconstructing patchy tau from the EB quadratic estimator
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
!*    :tlm[lx,ly] (dcmplx) : 2D Fourier modes of patchy tau, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: E, B
  double complex, intent(out), dimension(nx,ny) :: tlm
  !internal
  integer :: nn(2)
  double precision, dimension(nx,ny) :: els, lmask
  double complex, dimension(nx,ny) :: wE, wB, alm, ei2p

  nn = (/nx,ny/)
  call elarrays_2d(nn,D,els=els,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  wE = lmask*fC*E*ei2p
  wB = lmask*B*conjg(ei2p)

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wE,nn,D,-1)
  alm = aimag(wB*wE)
  call dft(alm,nn,D,1)

  tlm = alm

end subroutine qeb


end module rec_tau


