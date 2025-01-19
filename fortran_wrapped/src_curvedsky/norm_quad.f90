!////////////////////////////////////////////////////!
! * Normalization of quadratic tau/rot reconstruction
!////////////////////////////////////////////////////!

module norm_quad
  use alkernel, only: kernels_lens, kernels_tau, kernels_rot, kernels_lenstau, get_lfac
  !from basic
  use delensing
  implicit none

  private kernels_lens, kernels_tau, kernels_rot, kernels_lenstau, get_lfac

contains


subroutine qtt(est,lmax,rlmin,rlmax,TT,OCT,Al,lfac)
!*  Normalization of reconstructed fields from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :TT [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, TT, OCT
  !f2py intent(out) Al
  !f2py depend(rlmax) TT, OCT
  !f2py depend(lmax) Al
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: TT, OCT
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: rL(2), l
  double precision, dimension(lmax) :: lk2
  double precision, dimension(3,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG
  !opt4py :: lfac = ''

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtt): observed cltt is zero'
  end do

  W(1,:) = 1d0 / OCT(rlmin:rlmax)
  W(2,:) = TT(rlmin:rlmax)**2 / OCT(rlmin:rlmax)
  W(3,:) = TT(rlmin:rlmax) / OCT(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call Kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'S0')
    call Kernels_lens(rL,W(3,:),W(3,:),SG(2,:,:),'G0')
  case('amp')
    call Kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'G0')
  case('src')
    call Kernels_tau(rL,W(1,:),W(1,:),SG(1,1,:),'S0')
    call Kernels_tau(rL,W(1,:),W(1,:),SG(2,1,:),'G0')
    SG = SG/4d0
  end select

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  Al(2,1) = 0d0

end subroutine qtt

subroutine qte(est,lmax,rlmin,rlmax,TE,OCT,OCE,Al,lfac)
!*  Normalization of reconstructed fields from the TE quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, TE, OCT, OCE
  !f2py intent(out) Al
  !f2py depend(rlmax) TE, OCT, OCE
  !f2py depend(lmax) Al
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: TE, OCT, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(6,rlmin:rlmax) :: W
  double precision, dimension(3,2,lmax) :: SG
  !opt4py :: lfac = ''

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qte): observed clee is zero'
  end do

  W(1,:) = 1d0/OCT(rlmin:rlmax)
  W(2,:) = TE(rlmin:rlmax)**2/OCE(rlmin:rlmax)
  
  W(3,:) = TE(rlmin:rlmax)/OCT(rlmin:rlmax)
  W(4,:) = TE(rlmin:rlmax)/OCE(rlmin:rlmax)
  
  W(5,:) = 1d0/OCE(rlmin:rlmax)
  W(6,:) = TE(rlmin:rlmax)**2/OCT(rlmin:rlmax)

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
  case('rot')
    call kernels_rot(rL,W(5,:),W(6,:),SG(3,1,:),'Sp')
  case('src')
    call kernels_tau(rL,W(1,:),W(5,:),SG(1,1,:),'S0')
    call kernels_tau(rL,W(1,:),W(5,:),SG(2,1,:),'Gc')
    call kernels_tau(rL,W(1,:),W(5,:),SG(3,1,:),'Sp')
    SG = SG/4d0
  end select
  SG(2,:,:) = 2d0*SG(2,:,:)
  
  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(2,1) = 0d0
  case('rot')
    Al(1,1) = 0d0
  end select

end subroutine qte

subroutine qtb(est,lmax,rlmin,rlmax,TE,OCT,OCB,Al,lfac)
!*  Normalization of reconstructed fields from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :TE [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalizations (1 is dummy except lens = 0 and curl = 1), with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, TE, OCT, OCB
  !f2py intent(out) Al
  !f2py depend(rlmax) TE, OCT, OCB
  !f2py depend(lmax) Al
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: TE, OCT, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(2,rlmin:rlmax) :: W
  double precision, dimension(2,lmax) :: SG
  !opt4py :: lfac = ''

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtb): observed cltt is zero'
    if (OCB(l)==0d0) stop 'error (qtb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCB(rlmin:rlmax)
  W(2,:) = TE(rlmin:rlmax)**2 / OCT(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG,'Sm')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,:),'Sm')
  case('rot')
    call kernels_rot(rL,W(1,:),W(2,:),SG(1,:),'Sm')
  case('src')
    call kernels_tau(rL,W(1,:),W(1,:),SG(1,:),'Sm')
    SG = SG/4d0
  end select

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (SG(1,l)/=0d0)  Al(1,l) = lk2(l)/SG(1,l)
    if (SG(2,l)/=0d0)  Al(2,l) = lk2(l)/SG(2,l)
  end do
  select case(est)
  case('lens','amp','src')
    Al(1,1) = 0d0
  end select

end subroutine qtb

subroutine qee(est,lmax,rlmin,rlmax,EE,OCE,Al,lfac)
!*  Normalization of reconstructed amplitude modulation from the EE quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)   : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, EE, OCE
  !f2py intent(out) Al
  !f2py depend(rlmax) EE, OCE
  !f2py depend(lmax) Al
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE, OCE
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(3,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG
  !opt4py :: lfac = ''

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (qee): observed clee is zero'
  end do

  W(1,:) = 1d0/OCE(rlmin:rlmax)
  W(2,:) = EE(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  W(3,:) = EE(rlmin:rlmax) / OCE(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sp')
    call kernels_lens(rL,W(3,:),W(3,:),SG(2,:,:),'Gp')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('rot')
    call kernels_rot(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_rot(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('src')
    call kernels_tau(rL,W(1,:),W(1,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(1,:),W(1,:),SG(2,1,:),'Gp')
    SG = SG/4d0
  end select

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(2,1) = 0d0
  case('rot')
    Al(1,1) = 0d0
  end select

end subroutine qee

subroutine qeb(est,lmax,rlmin,rlmax,EE,OCE,OCB,BB,Al,lfac)
!*  Normalization of reconstructed fields from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :EE [l] (double)   : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double)  : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double)  : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optionals): 
!*    :BB [l] (double)   : Theory BB spectrum, with bounds (0:rlmax)
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, EE, BB, OCE, OCB
  !f2py intent(out) Al
  !f2py depend(rlmax) EE, BB, OCE, OCB
  !f2py depend(lmax) Al
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE, BB, OCE, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(6,rlmin:rlmax) :: W
  double precision, dimension(3,2,lmax) :: SG
  !opt4py :: lfac = ''
  !opt4py :: BB = 0
  !add2py :: if BB==0: BB=0.*EE

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (qeb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCE(rlmin:rlmax)
  W(2,:) = BB(rlmin:rlmax)**2 / OCB(rlmin:rlmax)
  W(3,:) = EE(rlmin:rlmax)/OCE(rlmin:rlmax)
  W(4,:) = BB(rlmin:rlmax)/OCB(rlmin:rlmax)
  W(5,:) = 1d0/OCB(rlmin:rlmax)
  W(6,:) = EE(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  
  SG = 0d0
  select case(est)
  case('lens')
    if (sum(BB)/=0d0) then
      call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sm')
      call kernels_lens(rL,W(3,:),W(4,:),SG(2,:,:),'Gm')
    end if
    call kernels_lens(rL,W(5,:),W(6,:),SG(3,:,:),'Sm')
  case('amp')
    if (sum(BB)/=0d0) then
      call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sm')
      call kernels_tau(rL,W(3,:),W(4,:),SG(2,1,:),'Gm')
    end if
    call kernels_tau(rL,W(5,:),W(6,:),SG(3,1,:),'Sm')
  case('rot')
    if (sum(BB)/=0d0) then
      call kernels_rot(rL,W(1,:),W(2,:),SG(1,1,:),'Sm')
      call kernels_rot(rL,W(3,:),W(4,:),SG(2,1,:),'Gm')
    end if
    call kernels_rot(rL,W(5,:),W(6,:),SG(3,1,:),'Sm')
  case('src')
    call kernels_tau(rL,W(1,:),W(5,:),SG(1,1,:),'Sm')
    call kernels_tau(rL,W(1,:),W(5,:),SG(2,1,:),'Gm')
    call kernels_tau(rL,W(1,:),W(5,:),SG(3,1,:),'Sm')
    SG = SG/4d0
  end select

  SG(2,:,:) = 2d0*SG(2,:,:)

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(1,1) = 0d0
  end select

end subroutine qeb

subroutine qbb(est,lmax,rlmin,rlmax,BB,OCB,Al,lfac)
!*  Normalization of reconstructed amplitude modulation from the BB quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :BB [l] (double)   : Theory BB spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double)  : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Al [2,l] (double) : Normalization, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, BB, OCB
  !f2py intent(out) Al
  !f2py depend(rlmax) BB, OCB
  !f2py depend(lmax) Al
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: BB, OCB
  double precision, intent(out), dimension(2,0:lmax) :: Al
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(3,rlmin:rlmax) :: W
  double precision, dimension(2,2,lmax) :: SG
  !opt4py :: lfac = ''

  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCB(l)==0d0) stop 'error (qbb): observed clbb is zero'
  end do

  W(1,:) = 1d0/OCB(rlmin:rlmax)
  W(2,:) = BB(rlmin:rlmax)**2 / OCB(rlmin:rlmax)
  W(3,:) = BB(rlmin:rlmax) / OCB(rlmin:rlmax)

  SG = 0d0
  select case(est)
  case('lens')
    call kernels_lens(rL,W(1,:),W(2,:),SG(1,:,:),'Sp')
    call kernels_lens(rL,W(3,:),W(3,:),SG(2,:,:),'Gp')
  case('amp')
    call kernels_tau(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('rot')
    call kernels_rot(rL,W(1,:),W(2,:),SG(1,1,:),'Sp')
    call kernels_rot(rL,W(3,:),W(3,:),SG(2,1,:),'Gp')
  case('src')
    call kernels_tau(rL,W(1,:),W(1,:),SG(1,1,:),'Sp')
    call kernels_tau(rL,W(1,:),W(1,:),SG(2,1,:),'Gp')
    SG = SG/4d0
  end select

  call get_lfac(lmax,lfac,lk2)

  Al = 0d0
  do l = 1, lmax
    if (sum(SG(:,1,l))/=0d0)  Al(1,l) = lk2(l)/sum(SG(:,1,l))
    if (sum(SG(:,2,l))/=0d0)  Al(2,l) = lk2(l)/sum(SG(:,2,l))
  end do
  select case(est)
  case('lens','amp','src')
    Al(2,1) = 0d0
  case('rot')
    Al(1,1) = 0d0
  end select

end subroutine qbb

subroutine qttte(est,lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TT and TE quadratic estimators
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCTT [l] (double): Theory TT spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double): Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, fCTT, fCTE, OCT, OCE, OCTE
  !f2py intent(out) Ig, Ic
  !f2py depend(rlmax) fCTT, fCTE, OCT, OCE, OCTE
  !f2py depend(lmax) Ig, Ic
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT, fCTE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3, W4, W5, W6, W7, W8
  double precision, dimension(2,lmax) :: S0, Gc, G0, Sc
  !opt4py :: lfac = ''

  write(*,*) 'norm qTTTE'

  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_qttte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_qttte): observed clee is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = 1d0/OCT(l)
    W2(l) = fCTT(l)*fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W3(l) = fCTE(l)/OCT(l)
    W4(l) = fCTT(l)*OCTE(l)/(OCT(l)*OCE(l))
    W5(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W6(l) = fCTT(l)/OCT(l)
    W7(l) = OCTE(l)/(OCT(l)*OCE(l))
    W8(l) = fCTT(l)*fCTE(l)/OCT(l)
  end do

  select case(est)
  case('lens')
    call kernels_lens(rL,W1,W2,S0,'S0')
    call kernels_lens(rL,W3,W4,Gc,'Gc')
    call kernels_lens(rL,W5,W6,G0,'G0')
    call kernels_lens(rL,W7,W8,Sc,'Sc')
  end select

  call get_lfac(lmax,lfac,lk2)

  do l = 1, lmax
    Ig(l) = S0(1,l)+Gc(1,l)+G0(1,l)+Sc(1,l)
    Ic(l) = S0(2,l)+Gc(2,l)+G0(2,l)+Sc(2,l)
    Ig(l) = Ig(l)/lk2(l)
    Ic(l) = Ic(l)/lk2(l)
  end do

end subroutine qttte

subroutine qttee(est,lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TT and EE quadratic estimators
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCTT [l] (double): Theory TT spectrum, with bounds (0:rlmax)
!*    :fCEE [l] (double): Theory EE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, fCTT, fCEE, OCT, OCE, OCTE
  !f2py intent(out) Ig, Ic
  !f2py depend(rlmax) fCTT, fCEE, OCT, OCE, OCTE
  !f2py depend(lmax) Ig, Ic
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT, fCEE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3, W4
  double precision, dimension(2,lmax) :: Sc, Gc
  !opt4py :: lfac = ''

  write(*,*) 'norm qTTEE'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0
  W1 = 0d0
  W2 = 0d0

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qttee): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qttee): observed clee is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W3(l) = fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W4(l) = fCTT(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do

  select case(est)
  case('lens')
    call kernels_lens(rL,W1,W2,Sc,'Sc')
    call kernels_lens(rL,W3,W4,Gc,'Gc')
  end select

  call get_lfac(lmax,lfac,lk2)

  do l = 1, lmax
    Ig(l) = Sc(1,l) + Gc(1,l)
    Ig(l) = Ig(l)/lk2(l)
  end do

end subroutine qttee

subroutine qteee(est,lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TE and EE quadratic estimators
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCEE [l] (double): Theory EE spectrum, with bounds (0:rlmax)
!*    :fCTE [l] (double): Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, fCEE, fCTE, OCT, OCE, OCTE
  !f2py intent(out) Ig, Ic
  !f2py depend(rlmax) fCEE, fCTE, OCT, OCE, OCTE
  !f2py depend(lmax) Ig, Ic
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCTE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3, W4, W5, W6, W7, W8
  double precision, dimension(2,lmax) :: Sc,Gp,Gc,Sp
  !opt4py :: lfac = ''

  write(*,*) 'norm qTEEE'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qteee): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qteee): observed clee is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTE(l)*fCEE(l)/OCE(l)
    W3(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W4(l) = fCEE(l)/OCE(l)
    W5(l) = fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W6(l) = fCTE(l)/OCE(l)
    W7(l) = 1d0/OCE(l)
    W8(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do

  select case(est)
  case('lens')
    call kernels_lens(rL,W1,W2,Sc,'Sc')
    call kernels_lens(rL,W3,W4,Gp,'Gp')
    call kernels_lens(rL,W5,W6,Gc,'Gc')
    call kernels_lens(rL,W7,W8,Sp,'Sp')
  end select

  call get_lfac(lmax,lfac,lk2)

  do l = 1, lmax
    Ig(l) = (Sc(1,l)+Gp(1,l)+Gc(1,l)+Sp(1,l))
    Ic(l) = (Sc(2,l)+Gp(2,l)+Gc(2,l)+Sp(2,l))
    Ig(l) = Ig(l)/lk2(l)
    Ic(l) = Ic(l)/lk2(l)
  end do

end subroutine qteee

subroutine qtbeb(est,lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,lfac)
!*  Correlation between unnormalized TB and EB quadratic estimators
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fCEE [l] (double): Theory EE spectrum, with bounds (0:rlmax)
!*    :fCBB [l] (double): Theory BB spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*    :OCTE [l] (double): Observed TE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, fCEE, fCBB, fCTE, OCT, OCE, OCTE, OCB
  !f2py intent(out) Ig, Ic
  !f2py depend(rlmax) fCEE, fCBB, fCTE, OCT, OCE, OCTE, OCB
  !f2py depend(lmax) Ig, Ic
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCBB,fCTE,OCT,OCE,OCTE,OCB
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(2,lmax) :: Gm, Sm
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3, W4
  !opt4py :: lfac = ''

  write(*,*) 'norm qTBEB'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (qtbeb): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (qtbeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (qtbeb): observed clbb is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCBB(l)/OCB(l)
    W3(l) = 1d0/OCB(l)
    W4(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do

  select case(est)
  case('lens')
    call kernels_lens(rL,W1,W2,Gm,'Gm')
    call kernels_lens(rL,W3,W4,Sm,'Sm')
  end select

  call get_lfac(lmax,lfac,lk2)

  do l = 1, lmax
    Ig(l) = Sm(1,l) + Gm(1,l)
    Ic(l) = Sm(2,l) + Gm(2,l)
    Ig(l) = Ig(l)/lk2(l)
    Ic(l) = Ic(l)/lk2(l)
  end do

end subroutine qtbeb

subroutine qmv(lmax,QDO,Al,Il,MV,Nl)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :lmax (int):    Maximum multipole of the output power spectra
!*    :QDO[6] (bool): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. Currently, BB is always False
!*    :Al [5,l] (double): Normalizations of each estimator (TT, TE, EE, TB, EB). 
!*    :Il [4,l] (double): Correlation between different estimators (TTxTE, TTxEE, TExEE, TBxEB).
!*
!*  Returns:
!*    :MV [l] (double):   Normalization of the MV estimator, with bounds (0:lmax)
!*    :Nl [6,l] (double): Weights for each estimator (TT, TE, EE, TB, EB, BB=0), with bounds (0:lmax)
!*
  !f2py intent(in) lmax, QDO, Al, Il
  !f2py intent(out) MV, Nl
  !f2py depend(lmax) Al, Il, MV, Nl
  implicit none
  ![input]
  integer, intent(in) :: lmax
  logical, intent(in), dimension(6) :: QDO
  double precision, intent(in), dimension(5,0:lmax) :: Al
  double precision, intent(in), dimension(4,0:lmax) :: Il
  double precision, intent(out), dimension(0:lmax) :: MV
  double precision, intent(out), dimension(6,0:lmax) :: Nl
  !internal
  integer :: qn = 5, QTT = 1, QTE = 2, QTB = 4, QEE = 3, QEB = 5!, QBB = 6
  integer :: X, Y, qmax, i, id(6), l
  double precision :: d1, d2, M1(5,5)
  double precision, allocatable :: M(:,:)

  write(*,*) 'norm qMV'
  id = 0

  qmax = qn
 
  MV = 0d0
  Nl = 0d0

  do l = 2, lmax

    ! noise covariance
    allocate(M(qmax,qmax));  M = 0d0

    if (QDO(QTT).and.QDO(QTE)) M(1,2) = Il(1,l)*Al(QTT,l)*Al(QTE,l)
    if (QDO(QTT).and.QDO(QEE)) M(1,3) = Il(2,l)*Al(QTT,l)*Al(QEE,l)
    if (QDO(QTE).and.QDO(QEE)) M(2,3) = Il(3,l)*Al(QTE,l)*Al(QEE,l)
    if (QDO(QTB).and.QDO(QEB)) M(4,5) = Il(4,l)*Al(QTB,l)*Al(QEB,l)
    do X = 1, qn
      M(X,X) = 1d10 ! some large value for not used estimator
      if (QDO(X)) M(X,X) = Al(X,l)
      do Y = X + 1, qn
        if(QDO(X).and.QDO(Y)) M(Y,X) = M(X,Y)
      end do
    end do

    ! inverting M with explict expression (symmetric)
    d1 = M(1,1)*(M(2,2)*M(3,3)-M(2,3)**2) + M(1,2)*(M(2,3)*M(3,1)-M(2,1)*M(3,3)) + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
    d2 = M(4,4)*M(5,5) - M(4,5)**2
    M1 = 0d0
    ! components
    M1(1,1) = M(2,2)*M(3,3) - M(2,3)**2
    M1(1,2) = M(1,3)*M(3,2) - M(1,2)*M(3,3)
    M1(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
    M1(2,2) = M(1,1)*M(3,3) - M(1,3)**2
    M1(2,3) = M(1,3)*M(2,1) - M(1,1)*M(2,3)
    M1(3,3) = M(1,1)*M(2,2) - M(1,2)**2
    M1(4,4) = M(5,5)
    M1(4,5) = -M(4,5)
    M1(5,5) = M(4,4)
    ! symmetric
    M1(2,1) = M1(1,2)
    M1(3,1) = M1(1,3)
    M1(3,2) = M1(2,3)
    M1(5,4) = M1(4,5)
    ! final
    M(1:3,1:3) = M1(1:3,1:3)/d1
    M(4:5,4:5) = M1(4:5,4:5)/d2

    MV(l)   = 1d0/sum(M)
    Nl(:,l) = sum(M,dim=2)

    deallocate(M)

  end do

end subroutine qmv

subroutine qall(est,QDO,lmax,rlmin,rlmax,fC,OC,Ag,Ac,Nlg,Nlc,lfac)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :QDO[6] (bool): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. 
!*    :lmax (int):    Maximum multipole of the output power spectra
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC/OC [l] (double): Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax) 
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ag [6,l] (double)  : Normalization of the TT, TE, EE, TB, EB, and MV estimators for lensing potential, with bounds (6,0:lmax)
!*    :Ac [6,l] (double)  : Same as Ag but for curl mode
!*    :Nlg [6,l] (double) : Weights for TT, TE, EE, TB, EB, and BB (=0) estimators for lensing potential, with bounds (6,0:lmax)
!*    :Nlc [6,l] (double) : Same as Nlg but for curl mode
!*
  !f2py intent(in) est, lfac, QDO, rlmin, rlmax, lmax, fC, OC
  !f2py intent(out) Ag, Ac, Nlg, Nlc
  !f2py depend(rlmax) fC, OC
  !f2py depend(lmax) Ag, Ac, Nlg, Nlc
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  logical, intent(in), dimension(6) :: QDO
  integer, intent(in) :: rlmin, rlmax, lmax
  double precision, intent(in), dimension(4,0:rlmax) :: fC, OC
  double precision, intent(out), dimension(6,0:lmax) :: Ag, Ac, Nlg, Nlc
  !internal
  character(1) :: gt
  integer :: TT = 1, EE = 2, BB = 3, TE = 4
  double precision, dimension(:,:), allocatable :: Ilg, Ilc
  double precision, dimension(6,2,0:lmax) :: Al
  !opt4py :: lfac = ''

  !//// interface ////!
  Ag  = 0d0
  Ac  = 0d0
  Nlg = 0d0
  Nlc = 0d0
  if (QDO(1))  call qtt(est,lmax,rlmin,rlmax,fC(TT,:),OC(TT,:),Al(1,:,:),lfac)
  if (QDO(2))  call qte(est,lmax,rlmin,rlmax,fC(TE,:),OC(TT,:),OC(EE,:),Al(2,:,:),lfac)
  if (QDO(3))  call qee(est,lmax,rlmin,rlmax,fC(EE,:),OC(EE,:),Al(3,:,:),lfac)
  if (QDO(4))  call qtb(est,lmax,rlmin,rlmax,fC(TE,:),OC(TT,:),OC(BB,:),Al(4,:,:),lfac)
  if (QDO(5))  call qeb(est,lmax,rlmin,rlmax,fC(EE,:),OC(EE,:),OC(BB,:),fC(BB,:),Al(5,:,:),lfac)
  Ag = Al(:,1,:)
  Ac = Al(:,2,:)

  allocate(Ilg(4,0:lmax),Ilc(4,0:lmax))
  if (QDO(1).and.QDO(2))  call qttte(est,lmax,rlmin,rlmax,fC(TT,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(1,:),Ilc(1,:),lfac)
  if (QDO(1).and.QDO(3))  call qttee(est,lmax,rlmin,rlmax,fC(TT,:),fC(EE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(2,:),Ilc(2,:),lfac)
  if (QDO(2).and.QDO(3))  call qteee(est,lmax,rlmin,rlmax,fC(EE,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(3,:),Ilc(3,:),lfac)
  if (QDO(4).and.QDO(5))  call qtbeb(est,lmax,rlmin,rlmax,fC(EE,:),fC(BB,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(BB,:),OC(TE,:),Ilg(4,:),Ilc(4,:),lfac)
  call qmv(lmax,QDO,Ag(1:5,0:lmax),Ilg,Ag(6,0:lmax),Nlg)
  call qmv(lmax,QDO,Ac(1:5,0:lmax),Ilc,Ac(6,0:lmax),Nlc)
  deallocate(Ilg,Ilc)

end subroutine qall

subroutine qeb_iter(lmax,elmax,rlmin,rlmax,dlmin,dlmax,CE,OCE,OCB,Cpp,Ag,Ac,iter,conv)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization
!*    :elmax (int)      : Maximum multipole of input EE spectra, CE and OCE
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :dlmin/dlmax (int): Minimum/Maximum multipole of E mode and lensing potential for delensing
!*    :CE [l] (double)  : Theory EE angular power spectrum, with bounds (0:elmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:elmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*    :Cpp [l] (double) : Theory lensing potential spectrum, with bounds (0:dlmax)
!*
!*  Args(optional):
!*    :iter (int)    : number of iteration, default to 1 (no iteration)
!*    :conv (double) : a parameter for convergence the iteration, default to 0.00001
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  !f2py intent(in) lmax, elmax, rlmin, rlmax, dlmin, dlmax, CE, OCE, OCB, Cpp, iter, conv
  !f2py intent(out) Ag, Ac
  !f2py depend(elmax) CE, OCE
  !f2py depend(rlmax) OCB
  !f2py depend(dlmax) Cpp
  !f2py depend(lmax) Ag, Ac
  implicit none
  !I/O
  integer, intent(in) :: lmax, elmax, rlmin, rlmax, dlmin, dlmax
  double precision, intent(in), dimension(0:elmax) :: CE, OCE
  double precision, intent(in), dimension(0:rlmax) :: OCB
  double precision, intent(in), dimension(0:dlmax) :: Cpp
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  integer, intent(in) :: iter
  double precision, intent(in) :: conv
  !internal
  integer :: i, n, l
  double precision :: ratio
  double precision :: Al(2,0:dlmax), rCBB(0:rlmax), BB(0:rlmax)
  !opt4py :: iter = 1
  !opt4py :: conv = 0.00001

  if (elmax<dlmax.or.elmax<rlmax) stop 'error (qeb_iter): does not support elmax<dlmax or elmax<rlmax'

  !initial values
  ratio = 1d0
  rCBB  = OCB
  Al   = 0d0
  BB   = 0d0

  do n = 1, iter !loop for iteration 

    !lensing reconstruction with EB
    call qeb('lens',dlmax,rlmin,rlmax,CE(0:rlmax),OCE(0:rlmax),rCBB,BB(0:rlmax),Al(:,0:dlmax),'')

    !convergence check using gradient mode
    if (n>=2) then
      ratio = (sum(Ag)/sum(Al(1,:))-1d0)/dble(dlmax)
      write(*,*) n, ratio
    end if
    Ag = Al(1,:)
    Ac = Al(2,:)

    if (abs(ratio) < conv) exit

    !delensing with EB-estimator
    call clbb_est((/rlmin,rlmax/),(/dlmin,dlmax/),(/dlmin,dlmax/),CE(1:dlmax),Cpp(1:dlmax),OCE(1:dlmax)-CE(1:dlmax),Al(1,1:dlmax),rCBB(1:rlmax))
    rCBB = OCB - rCBB !delensed B-mode

    if(n==iter) write(*,*) 'not well converged'

  end do

end subroutine qeb_iter

subroutine xtt(est,lmax,rlmin,rlmax,fC,OCT,Ag,lfac)
!*  Unnormalized response between lensing potential and amplitude modulation from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :lfac (str)       : Multiplying square of L(L+1)/2, i.e., convergence (lfac='k') or lensing potential (lfac='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : Cross normalization, with bounds (0:lmax)
!*
  !f2py intent(in) est, lfac, lmax, rlmin, rlmax, fC, OCT
  !f2py intent(out) Ag
  !f2py depend(rlmax) fC, OCT
  !f2py depend(lmax) Ag
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(1), intent(in) :: lfac
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3
  double precision, dimension(2,lmax) :: S0, G0
  !opt4py :: lfac = ''

  write(*,*) 'cross norm TT'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_xtt): observed cltt is zero'
  end do

  W1 = 1d0 / OCT(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2 / OCT(rlmin:rlmax)
  W3 = fC(rlmin:rlmax) / OCT(rlmin:rlmax)

  S0 = 0d0
  G0 = 0d0
  select case(est)
  case('lensamp')
    call kernels_lenstau(rL,W1,W2,S0,'S0')
    call kernels_lenstau(rL,W3,W3,G0,'G0')
  case('lenssrc')
    call kernels_lenstau(rL,W1,W3,S0,'S0')
    call kernels_lenstau(rL,W1,W3,G0,'G0')
    S0 = S0*0.5d0
    G0 = G0*0.5d0
  case('ampsrc')
    call Kernels_tau(rL,W1,W3,S0(1,:),'S0')
    call Kernels_tau(rL,W1,W3,G0(1,:),'G0')
    S0 = S0*0.5d0
    G0 = G0*0.5d0
  end select

  call get_lfac(lmax,lfac,lk2)

  Ag = 0d0
  do l = 1, lmax
    Ag(l) = (S0(1,l)+G0(1,l))/dsqrt(lk2(l))
  end do

end subroutine xtt

subroutine xeb(est,lmax,rlmin,rlmax,EE,EB,OCE,OCB,BB,Aa)
!*  Response of reconstructed field to other in the EB quadratic estimator
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
  !f2py intent(in) est, lmax, rlmin, rlmax, EE, EB, OCE, OCB, BB
  !f2py intent(out) Aa
  !f2py depend(rlmax) EE, EB, OCE, OCB, BB
  !f2py depend(lmax) Aa
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: EE, EB, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Aa
  !optional
  double precision, intent(in), dimension(0:rlmax), optional :: BB
  !f2py double precision :: BB = 0
  !docstr :: BB = EE*0
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2, W3, W4
  double precision, dimension(3,lmax) :: SG

  rL = (/rlmin,rlmax/)
  SG = 0d0

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_quad.xeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_quad.xeb): observed clbb is zero'
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
  W3 = 1d0/OCB(rlmin:rlmax)
  W4 = EE(rlmin:rlmax)*EB(rlmin:rlmax) / OCE(rlmin:rlmax)

  select case(est)
  case('rotamp')
    call kernels_rot(rL,W1,W2,SG(2,:),'Gm')
    call kernels_rot(rL,W3,W4,SG(3,:),'Sm')
  end select

  Aa = sum(SG,dim=1)/2d0

end subroutine xeb




end module norm_quad




