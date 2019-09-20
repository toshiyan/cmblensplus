!////////////////////////////////////////////////////!
! Normalization of quadratic rotation reconstruction
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
!*    :fC [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
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


subroutine qeb(lmax,rlmin,rlmax,EE,OCE,OCB,BB,Aa)
!*  Normalization of reconstructed pol. rot. angle from the EB quadratic estimator
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
!*    :Aa [l] (double) : Pol. rot. angle normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Aa
  !optional
  double precision, intent(in), dimension(0:rlmax), optional :: BB
  !f2py double precision :: BB = 0
  !docstr :: BB = EE*0
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(3,lmax) :: SG

  write(*,*) 'norm qEB (rot)'
  rL = (/rlmin,rlmax/)
  SG = 0d0

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_rot.qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_rot.qeb): observed clbb is zero'
  end do

  if (present(BB).and.sum(BB)/=0d0) then

    W1 = 1d0/OCE(rlmin:rlmax)
    W2 = BB(rlmin:rlmax)**2 / OCB(rlmin:rlmax)
    call kernels_rot(rL,W1,W2,SG(1,:),'Sm')

    W1 = EE(rlmin:rlmax)/OCE(rlmin:rlmax)
    W2 = BB(rlmin:rlmax)/OCB(rlmin:rlmax)
    call kernels_rot(rL,W1,W2,SG(2,:),'Gm')
    SG(2,:) = -2d0*SG(2,:)

  end if

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = EE(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  call kernels_rot(rL,W1,W2,SG(3,:),'Sm')

  Aa = 0d0
  do l = 1, lmax
    if (sum(SG(:,l))/=0d0)  Aa(l) = 1d0/sum(SG(:,l))
  end do

end subroutine qeb


subroutine teb(lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB,Aa)
!*  Response of reconstructed pol. rot. angle to amplitude in the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :EB [l] (double)  : Theory EB spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optionals): 
!*    :BB [l] (double)  : Theory BB spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :Aa [l] (double) : Pol. rot. angle normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE, EB, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Aa
  !optional
  double precision, intent(in), dimension(0:rlmax), optional :: BB
  !f2py double precision :: BB = 0
  !docstr :: BB = EE*0
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(3,lmax) :: SG

  write(*,*) 'norm qEB (rot)'
  rL = (/rlmin,rlmax/)
  SG = 0d0

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_rot.qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_rot.qeb): observed clbb is zero'
  end do

  if (present(BB).and.sum(BB)/=0d0) then

    W1 = 1d0/OCE(rlmin:rlmax)
    W2 = BB(rlmin:rlmax)*EB(rlmin:rlmax) / OCB(rlmin:rlmax)
    call kernels_rot(rL,W1,W2,SG(1,:),'Sm')
    SG(1,:) = -SG(1,:)

  end if

  W1 = EE(rlmin:rlmax)/OCE(rlmin:rlmax)
  W2 = EB(rlmin:rlmax)/OCB(rlmin:rlmax)
  !W2 = (EB(rlmin:rlmax)-BB(rlmin:rlmax))/OCB(rlmin:rlmax)
  call kernels_rot(rL,W1,W2,SG(2,:),'Gm')

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = EE(rlmin:rlmax)*EB(rlmin:rlmax) / OCE(rlmin:rlmax)
  call kernels_rot(rL,W1,W2,SG(3,:),'Sm')

  Aa = sum(SG,dim=1)/2d0

end subroutine teb



end module norm_rot


