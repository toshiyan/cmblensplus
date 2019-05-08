!////////////////////////////////////////////////////!
! * Normalization of quadratic tau reconstruction
!////////////////////////////////////////////////////!

module norm_tau
  use alkernel, only: Kernels_Tau
  implicit none

  private Kernels_Tau

contains


subroutine qtt(lmax,rlmin,rlmax,fC,OCT,At)
!*  Normalization of reconstructed tau from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory temperature angular power spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed temperature angular power spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :At [l] (double) : tau normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: At
  !internal
  integer :: rL(2), l
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(min(2*rlmax,lmax)) :: S0, G0

  write(*,*) 'norm qTT (tau)'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_tau.qtt): observed cltt is zero'
  end do

  !filtering functions
  W1 = 1d0 / OCT(rlmin:rlmax)

  !main calculation
  W2 = W1 * fC(rlmin:rlmax)**2
  S0 = 0d0
  call Kernels_tau(rL,W1,W2,S0,'S0')

  W2 = W1 * fC(rlmin:rlmax)
  G0 = 0d0
  call Kernels_tau(rL,W2,W2,G0,'G0')

  At = 0d0
  do l = 1, lmax
    if (S0(l)+G0(l)/=0d0)  At(l) = 1d0/(S0(l)+G0(l))
  end do

end subroutine qtt


end module norm_tau

