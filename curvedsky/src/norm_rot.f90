!////////////////////////////////////////////////////!
! * Normalization of quadratic rotation reconstruction
!////////////////////////////////////////////////////!

module norm_rot
  use alkernel,   only: kernels_rot
  implicit none

  private kernels_rot

contains


subroutine qeb(lmax,rlmin,rlmax,fCEE,OCEE,OCBB,Al)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE, OCEE, OCBB
  double precision, intent(out), dimension(0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(lmax) :: Sm

  write(*,*) 'norm qEB (rot)'

  rL = (/rlmin,rlmax/)
  W1 = 1d0/OCBB(rlmin:rlmax)
  W2 = fCEE(rlmin:rlmax)**2 / OCEE(rlmin:rlmax)
  Sm = 0d0
  call kernels_rot(rL,W1,W2,Sm,'Sm')

  Al = 0d0
  do l = 1, lmax
    if (Sm(l)/=0d0)  Al(l) = 1d0/Sm(l)
  end do

end subroutine qeb

end module norm_rot

