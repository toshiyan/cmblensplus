!////////////////////////////////////////////////////!
! * Normalization of quadratic tau reconstruction
!////////////////////////////////////////////////////!

module norm_tau
  use alkernel, only: kernels_tau, kernels_rot
  implicit none

  private kernels_tau, kernels_rot

contains


subroutine qtt(lmax,rlmin,rlmax,fC,OCT,At)
!*  Normalization of reconstructed amplitude modulation from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
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


subroutine qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB,At)
!*  Normalization of reconstructed amplitude modulation from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optionals): 
!*    :BB [l] (double)  : Theory BB spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :At [l] (double)  : Amplitude modulation normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: At
  !optional
  double precision, intent(in), dimension(0:rlmax), optional :: BB
  !f2py double precision :: BB = 0
  !docstr :: BB = EE*0
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(3,lmax) :: SG

  write(*,*) 'norm qEB (tau)'
  rL = (/rlmin,rlmax/)
  SG = 0d0

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_rot.qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_rot.qeb): observed clbb is zero'
  end do

  if (present(BB).and.sum(BB)/=0d0) then

    W1 = 1d0/OCE(rlmin:rlmax)
    W2 = BB(rlmin:rlmax)**2 / OCB(rlmin:rlmax)
    call kernels_tau(rL,W1,W2,SG(1,:),'Sm')

    W1 = EE(rlmin:rlmax)/OCE(rlmin:rlmax)
    W2 = BB(rlmin:rlmax)/OCB(rlmin:rlmax)
    call kernels_tau(rL,W1,W2,SG(2,:),'Gm')
    SG(2,:) = 2d0*SG(2,:)

  end if

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = EE(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  call kernels_tau(rL,W1,W2,SG(3,:),'Sm')

  At = 0d0
  do l = 1, lmax
    if (sum(SG(:,l))/=0d0)  At(l) = 1d0/sum(SG(:,l))
  end do

end subroutine qeb


subroutine oeb(lmax,rlmin,rlmax,EB,OCE,OCB,At)
!*  Normalization of reconstructed amplitude from the EB quadratic estimator
!*  The kernels are the same as the rotation normalization but with a factor 4. 
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :EB [l] (double)  : Theory EB spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :At [l] (double)  : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EB, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: At
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(3,lmax) :: SG

  write(*,*) 'norm qEB (tau)'
  rL = (/rlmin,rlmax/)
  SG = 0d0

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_rot.qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_rot.qeb): observed clbb is zero'
  end do

  W1 = 1d0/OCE(rlmin:rlmax)
  W2 = EB(rlmin:rlmax)**2 / OCB(rlmin:rlmax)
  call kernels_rot(rL,W1,W2,SG(1,:),'Sm')

  W1 = EB(rlmin:rlmax)/OCE(rlmin:rlmax)
  W2 = EB(rlmin:rlmax)/OCB(rlmin:rlmax)
  call kernels_rot(rL,W1,W2,SG(2,:),'Gm')
  SG(2,:) = 2d0*SG(2,:)

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = EB(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  call kernels_rot(rL,W1,W2,SG(3,:),'Sm')

  At = 0d0
  do l = 1, lmax
    if (sum(SG(:,l))/=0d0)  At(l) = 4d0/sum(SG(:,l))
  end do

end subroutine oeb



end module norm_tau

