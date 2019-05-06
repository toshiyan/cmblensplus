!////////////////////////////////////////////////////!
! * Normalization of quadratic tau reconstruction
!////////////////////////////////////////////////////!

module norm_tau
  use alkernel, only: Kernels_Tau
  implicit none

  private Kernels_Tau

contains


subroutine qtt(lmax,rlmin,rlmax,fC,OCTT,Al)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCTT
  double precision, intent(out), dimension(0:lmax) :: Al
  !internal
  integer :: oL(2), rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(min(2*rlmax,lmax)) :: S0, G0

  write(*,*) 'norm qTT (tau)'

  rL = (/rlmin,rlmax/)

  !only compute nonzero values
  oL = [1,min(2*rlmax,lmax)]

  !filtering functions
  W1 = 1d0 / OCTT(rlmin:rlmax)

  !main calculation
  W2 = W1 * fC(rlmin:rlmax)**2
  S0 = 0d0
  call Kernels_tau(rL,W1,W2,S0,'S0')

  W2 = W1 * fC(rlmin:rlmax)
  G0 = 0d0
  call Kernels_tau(rL,W2,W2,G0,'G0')

  Al = 0d0
  Al(oL(1):oL(2)) = 1d0/(S0(oL(1):oL(2))+G0(oL(1):oL(2)))

end subroutine qtt

end module norm_tau

