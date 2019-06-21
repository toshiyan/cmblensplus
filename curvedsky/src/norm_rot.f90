!////////////////////////////////////////////////////!
! * Normalization of quadratic rotation reconstruction
!////////////////////////////////////////////////////!

module norm_rot
  use alkernel,   only: kernels_rot
  implicit none

  private kernels_rot

contains


subroutine qtb(lmax,rlmin,rlmax,fC,OCT,OCB,Aa)
!*  Normalization of reconstructed pol. rot. angle from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory TE angular power spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed temperature angular power spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed B-mode angular power spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :Aa [l] (double) : Pol. rot. angle normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT, OCB
  double precision, intent(out), dimension(0:lmax) :: Aa
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(lmax) :: Sm

  write(*,*) 'norm qTB (rot)'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_rot.qtb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_rot.qtb): observed clbb is zero'
  end do

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2 / OCT(rlmin:rlmax)
  Sm = 0d0
  call kernels_rot(rL,W1,W2,Sm,'Sm')

  Aa = 0d0
  do l = 1, lmax
    if (Sm(l)/=0d0)  Aa(l) = 1d0/Sm(l)
  end do

end subroutine qtb


subroutine qeb(lmax,rlmin,rlmax,fC,OCE,OCB,Aa)
!*  Normalization of reconstructed pol. rot. angle from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory EE angular power spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed E-mode angular power spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed B-mode angular power spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :Aa [l] (double) : Pol. rot. angle normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Aa
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(lmax) :: Sm

  write(*,*) 'norm qEB (rot)'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_rot.qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_rot.qeb): observed clbb is zero'
  end do

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  Sm = 0d0
  call kernels_rot(rL,W1,W2,Sm,'Sm')

  Aa = 0d0
  do l = 1, lmax
    if (Sm(l)/=0d0)  Aa(l) = 1d0/Sm(l)
  end do

end subroutine qeb

end module norm_rot

