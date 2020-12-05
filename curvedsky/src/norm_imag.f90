!//////////////////////////////////////////////////////////!
! Normalization of quadratic imaginary lens reconstruction
!//////////////////////////////////////////////////////////!

module norm_imag
  use alkernel, only: get_lfac, kernels_lens, kernels_tau, kernels_lenstau
  implicit none

  private get_lfac, kernels_lens, kernels_tau, kernels_lenstau

contains


subroutine qte(est,lmax,rlmin,rlmax,TB,OCT,OCE,Al,lfac)
!*  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the TE quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TB [l] (double)  : Theory TB spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)    : Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: TB, OCT, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(2,rlmin:rlmax) :: W
  double precision, dimension(2,lmax) :: SG
  !opt4py :: lfac = ''

  write(*,*) 'norm qTE'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qte): observed clee is zero'
  end do

  W(1,:) = 1d0 / OCE(rlmin:rlmax)
  W(2,:) = TB(rlmin:rlmax)**2 / OCT(rlmin:rlmax)
  
  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG,'Sm')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,:),'Sm')
  end select

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (SG(1,l)/=0d0)  Al(1,l) = lk2(l)/SG(1,l)/4d0
    if (SG(2,l)/=0d0)  Al(2,l) = lk2(l)/SG(2,l)/4d0
  end do
  Al(1,1) = 0d0

end subroutine qte


subroutine qtb(est,lmax,rlmin,rlmax,TB,OCT,OCB,Al,lfac)
!*  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TB [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)    : Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: TB, OCT, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(6,rlmin:rlmax) :: W
  double precision, dimension(3,2,lmax) :: SG
  !opt4py :: lfac = ''

  write(*,*) 'norm qTB'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtb): observed cltt is zero'
    if (OCB(l)==0d0) stop 'error (qtb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCT(rlmin:rlmax)
  W(2,:) = TB(rlmin:rlmax)**2/OCB(rlmin:rlmax)
  W(3,:) = TB(rlmin:rlmax)/OCT(rlmin:rlmax)
  W(4,:) = TB(rlmin:rlmax)/OCB(rlmin:rlmax)
  W(5,:) = 1d0/OCB(rlmin:rlmax)
  W(6,:) = TB(rlmin:rlmax)**2/OCT(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'S0')
    call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gc')
    call kernels_lens(rL,W(5,:),W(6,:),SG(3,:,:),'Sp')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'S0')
    call kernels_tau(rL,W(3,:),W(4,:),SG(2,1,:),'Gc')
    call kernels_tau(rL,W(5,:),W(6,:),SG(3,1,:),'Sp')
  end select
  SG(2,:,:) = 2d0*SG(2,:,:)

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))/4d0
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))/4d0
  end do
  Al(2,1) = 0d0

end subroutine qtb


subroutine qee(est,lmax,rlmin,rlmax,fC,OCE,Al,lfac)
!*  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the E-mode quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)    : Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(3,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG
  !opt4py :: lfac = ''

  write(*,*) 'norm qEE'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (qee): observed clee is zero'
  end do

  W(1,:) = 1d0 / OCE(rlmin:rlmax)
  W(2,:) = fC(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  W(3,:) = fC(rlmin:rlmax) / OCE(rlmin:rlmax)
  
  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sm')
    call kernels_lens(rL,W(3,:),W(3,:),SG(2,:,:),'Gm')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sm')
    call kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'Gm')
  end select
  SG(2,:,:) = -SG(2,:,:)

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))/4d0
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))/4d0
  end do
  Al(1,1) = 0d0

end subroutine qee


subroutine qeb(est,lmax,rlmin,rlmax,fC,OCE,OCB,Al,lfac)
!*  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory EB spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)    : Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCE, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(6,rlmin:rlmax) :: W
  double precision, dimension(3,2,lmax) :: SG
  !opt4py :: lfac = ''

  write(*,*) 'norm qEB'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (qeb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCB(rlmin:rlmax)
  W(2,:) = fC(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  W(3,:) = fC(rlmin:rlmax)/OCE(rlmin:rlmax)
  W(4,:) = fC(rlmin:rlmax)/OCB(rlmin:rlmax)
  W(5,:) = 1d0/OCE(rlmin:rlmax)
  W(6,:) = fC(rlmin:rlmax)**2 / OCB(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sp')
    call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gp') !px Gamma
    call kernels_lens(rL,W(5,:),W(6,:),SG(3,:,:),'Sp')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(3,:),W(4,:),SG(2,1,:),'Gp') !px Gamma
    call kernels_tau(rL,W(5,:),W(6,:),SG(3,1,:),'Sp')
  end select
  SG(2,:,:) = 2d0*SG(2,:,:)

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))/4d0
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))/4d0
  end do
  Al(2,1) = 0d0

end subroutine qeb


subroutine qbb(est,lmax,rlmin,rlmax,fC,OCB,Al,lfac)
!*  Normalization of reconstructed imaginary CMB lensing potential and its curl mode from the B-mode quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory BB spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)    : Type of output, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(3,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG
  !opt4py :: lfac = ''

  write(*,*) 'norm qBB'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCB(l)==0d0) stop 'error (qbb): observed clbb is zero'
  end do

  W(1,:) = 1d0 / OCB(rlmin:rlmax)
  W(2,:) = fC(rlmin:rlmax)**2 / OCB(rlmin:rlmax)
  W(3,:) = fC(rlmin:rlmax) / OCB(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sm')
    call kernels_lens(rL,W(3,:),W(3,:),SG(2,:,:),'Gm')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sm')
    call kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'Gm')
  end select
  SG(2,:,:) = -SG(2,:,:)

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))/4d0
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))/4d0
  end do
  Al(1,1) = 0d0

end subroutine qbb


end module norm_imag


