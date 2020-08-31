!////////////////////////////////////////////////////!
! * Normalization of quadratic src reconstruction
!////////////////////////////////////////////////////!

module norm_src
  use alkernel, only: kernels_tau
  implicit none

  private kernels_tau

contains


subroutine qtt(lmax,rlmin,rlmax,OCT,As)
!*  Normalization of reconstructed src field from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :As [l] (double) : src field normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: OCT
  double precision, intent(out), dimension(0:lmax) :: As
  !internal
  integer :: rL(2), l
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(lmax) :: S0, G0

  write(*,*) 'norm qTT (src)'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_src.qtt): observed cltt is zero'
  end do

  W1 = 0.5d0 / OCT(rlmin:rlmax)
  S0 = 0d0
  call kernels_tau(rL,W1,W1,S0,'S0')

  G0 = 0d0
  call kernels_tau(rL,W1,W1,G0,'G0')

  As = 0d0
  do l = 1, lmax
    if (S0(l)+G0(l)/=0d0)  As(l) = 1d0/(S0(l)+G0(l))
  end do

end subroutine qtt


end module norm_src

