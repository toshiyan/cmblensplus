!////////////////////////////////////////////////////!
! * Normalization of quadratic src reconstruction
!////////////////////////////////////////////////////!

module rec_src
  use alkernel, only: kernels_tau
  implicit none

  private kernels_tau

contains


subroutine qtt(lmax,rlmin,rlmax,OCTT,Al)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: OCTT
  double precision, intent(out), dimension(0:lmax) :: Al
  !internal
  integer :: rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(lmax) :: S0, G0

  write(*,*) 'norm qTT (src)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0 / OCTT(rlmin:rlmax)
  S0 = 0d0
  call kernels_tau(rL,W1,W1,S0,'S0')

  G0 = 0d0
  call kernels_tau(rL,W1,W1,G0,'G0')

  Al(1:lmax) = 1d0/(S0(1:lmax)+G0(1:lmax))

end subroutine qtt


end module rec_src

