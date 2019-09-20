!////////////////////////////////////////////////!
! Normalization of quadratic lens reconstruction
!////////////////////////////////////////////////!

module norm_lens
  use alkernel, only: kernels_lens, kernels_lenstau
  use delensing
  !use lapack95, only: inv_lapack
  implicit none

  private kernels_lens, kernels_lenstau
  !private inv_lapack

contains


subroutine qtt(lmax,rlmin,rlmax,fC,OCT,Ag,Ac,gtype)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double)   : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, G0

  write(*,*) 'norm qTT (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qtt): observed cltt is zero'
  end do

  W1 = 1d0 / OCT(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  S0 = 0d0
  call kernels_lens(rL,W1,W2,S0,'S0')

  W2 = W1 * fC(rlmin:rlmax)
  G0 = 0d0
  call kernels_lens(rL,W2,W2,G0,'G0')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (S0(1,l)+G0(1,l)/=0d0)  Ag(l) = lk2(l)/(S0(1,l)+G0(1,l))
    if (S0(2,l)+G0(2,l)/=0d0)  Ac(l) = lk2(l)/(S0(2,l)+G0(2,l))
  end do
  Ac(1) = 0d0

end subroutine qtt


subroutine qte(lmax,rlmin,rlmax,fC,OCT,OCE,Ag,Ac,gtype)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the TE quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*

  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCT, OCE
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, Sp, Gc

  write(*,*) 'norm qTE (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_lens.qte): observed clee is zero'
  end do


  W1 = 1d0/OCT(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2/OCE(rlmin:rlmax)
  S0 = 0d0
  call kernels_lens(rL,W1,W2,S0,'S0')

  W1 = 1d0/OCE(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2/OCT(rlmin:rlmax)
  Sp = 0d0
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  W1 = fC(rlmin:rlmax)/OCT(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)/OCE(rlmin:rlmax)
  Gc = 0d0
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (S0(1,l)+Sp(1,l)+2d0*Gc(1,l)/=0d0)  Ag(l) = lk2(l)/(S0(1,l)+Sp(1,l)+2d0*Gc(1,l))
    if (S0(2,l)+Sp(2,l)+2d0*Gc(2,l)/=0d0)  Ac(l) = lk2(l)/(S0(2,l)+Sp(2,l)+2d0*Gc(2,l))
  end do
  Ac(1) = 0d0

end subroutine qte


subroutine qtb(lmax,rlmin,rlmax,fC,OCT,OCB,Ag,Ac,gtype)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the TB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory TE spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double) : Observed TT spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCT, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sm

  write(*,*) 'norm qTB (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qtb): observed cltt is zero'
    if (OCB(l)==0d0) stop 'error (norm_lens.qtb): observed clbb is zero'
  end do

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2/OCT(rlmin:rlmax)
  Sm = 0d0
  call kernels_lens(rL,W1,W2,Sm,'Sm')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (Sm(1,l)/=0d0)  Ag(l) = lk2(l)/Sm(1,l)
    if (Sm(2,l)/=0d0)  Ac(l) = lk2(l)/Sm(2,l)
  end do
  Ag(1) = 0d0

end subroutine qtb


subroutine qee(lmax,rlmin,rlmax,fC,OCE,Ag,Ac,gtype)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the E-mode quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCE
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sp, Gp

  write(*,*) 'norm qEE (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_lens.qee): observed clee is zero'
  end do

  W1 = 1d0/OCE(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  Sp = 0d0
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  W2 = W1 * fC(rlmin:rlmax)
  Gp = 0d0
  call kernels_lens(rL,W2,W2,Gp,'Gp')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (Sp(1,l)+Gp(1,l)/=0d0)  Ag(l) = lk2(l)/(Sp(1,l)+Gp(1,l))
    if (Sp(2,l)+Gp(2,l)/=0d0)  Ac(l) = lk2(l)/(Sp(2,l)+Gp(2,l))
  end do
  Ac(1) = 0d0

end subroutine qee


subroutine qeb(lmax,rlmin,rlmax,fC,OCE,OCB,Ag,Ac,gtype)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory EE spectrum, with bounds (0:rlmax)
!*    :OCE [l] (double) : Observed EE spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sm

  write(*,*) 'norm qEB (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCE(l)==0d0) stop 'error (norm_lens.qeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_lens.qeb): observed clbb is zero'
  end do

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  Sm = 0d0
  call kernels_lens(rL,W1,W2,Sm,'Sm')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (Sm(1,l)/=0d0)  Ag(l) = lk2(l)/Sm(1,l)
    if (Sm(2,l)/=0d0)  Ac(l) = lk2(l)/Sm(2,l)
  end do
  Ag(1) = 0d0

end subroutine qeb


subroutine qbb(lmax,rlmin,rlmax,fC,OCB,Ag,Ac,gtype)
!*  Normalization of reconstructed CMB lensing potential and its curl mode from the B-mode quadratic estimator
!*
!*  Args:
!*    :lmax (int)       : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int): Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)  : Theory BB spectrum, with bounds (0:rlmax)
!*    :OCB [l] (double) : Observed BB spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''

  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sp, Gp

  write(*,*) 'norm qBB (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCB(l)==0d0) stop 'error (norm_lens.qbb): observed clbb is zero'
  end do

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  W2 = W1 * fC(rlmin:rlmax)
  call kernels_lens(rL,W2,W2,Gp,'Gp')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if(Sp(1,l)+Gp(1,l)/=0d0) Ag(l) = lk2(l)/(Sp(1,l)+Gp(1,l))
    if(Sp(2,l)+Gp(2,l)/=0d0) Ac(l) = lk2(l)/(Sp(2,l)+Gp(2,l))
  end do
  Ac(1) = 0d0

end subroutine qbb


subroutine qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype)
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
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT, fCTE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, Gc, G0, Sc

  write(*,*) 'norm qTTTE (lens)'

  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qttte): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_lens.qttte): observed clee is zero'
  end do


  do l = rlmin, rlmax
    W1(l) = 1d0/OCT(l)
    W2(l) = fCTT(l)*fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,S0,'S0')

  do l = rlmin, rlmax
    W1(l) = fCTE(l)/OCT(l)
    W2(l) = fCTT(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)/OCT(l)
  end do
  call kernels_lens(rL,W1,W2,G0,'G0')

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*fCTE(l)/OCT(l)
  end do
  call kernels_lens(rL,W1,W2,Sc,'Sc')

  do l = 1, lmax
    Ig(l) = S0(1,l)+Gc(1,l)+G0(1,l)+Sc(1,l)
    Ic(l) = S0(2,l)+Gc(2,l)+G0(2,l)+Sc(2,l)
    Ig(l) = Ig(l)/lk2(l)
    Ic(l) = Ic(l)/lk2(l)
  end do

end subroutine qttte


subroutine qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic,gtype)
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
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT, fCEE, OCT, OCE, OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sc, Gc

  write(*,*) 'norm qTTEE (lens)'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0
  W1 = 0d0
  W2 = 0d0

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qttee): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_lens.qttee): observed clee is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Sc,'Sc')

  do l = rlmin, rlmax
    W1(l) = fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  do l = 1, lmax
    Ig(l) = Sc(1,l) + Gc(1,l)
    Ig(l) = Ig(l)/lk2(l)
  end do

end subroutine qttee


subroutine qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic,gtype)
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
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCTE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sc,Gp,Gc,Sp

  write(*,*) 'norm qTEEE (lens)'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qteee): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_lens.qteee): observed clee is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTE(l)*fCEE(l)/OCE(l)
  end do
  call kernels_lens(rL,W1,W2,Sc,'Sc')

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCEE(l)/OCE(l)
  end do
  call kernels_lens(rL,W1,W2,Gp,'Gp')

  do l = rlmin, rlmax
    W1(l) = fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTE(l)/OCE(l)
  end do
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  do l = rlmin, rlmax
    W1(l) = 1d0/OCE(l)
    W2(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  do l = 1, lmax
    Ig(l) = (Sc(1,l)+Gp(1,l)+Gc(1,l)+Sp(1,l))
    Ic(l) = (Sc(2,l)+Gp(2,l)+Gc(2,l)+Sp(2,l))
    Ig(l) = Ig(l)/lk2(l)
    Ic(l) = Ic(l)/lk2(l)
  end do

end subroutine qteee


subroutine qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic,gtype)
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
!*    :gtype (str)    : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ig [l] (double) : Correlation between lensing potential estimators, with bounds (0:lmax)
!*    :Ic [l] (double) : Correlation between curl mode estimators, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCBB,fCTE,OCT,OCE,OCTE,OCB
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(2,lmax) :: Gm, Sm
  double precision, dimension(rlmin:rlmax) :: W1, W2

  write(*,*) 'norm qTBEB (lens)'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.qtbeb): observed cltt is zero'
    if (OCE(l)==0d0) stop 'error (norm_lens.qtbeb): observed clee is zero'
    if (OCB(l)==0d0) stop 'error (norm_lens.qtbeb): observed clbb is zero'
  end do

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCBB(l)/OCB(l)
  end do
  call kernels_lens(rL,W1,W2,Gm,'Gm')

  do l = rlmin, rlmax
    W1(l) = 1d0/OCB(l)
    W2(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Sm,'Sm')

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

  write(*,*) 'norm qMV (lens)'
  id = 0

  qmax = qn
 
  MV = 0d0
  Nl = 0d0

  do l = 2, lmax

    ! noise covariance
    allocate(M(qmax,qmax));  M = 0d0

    !if (QDO(QTT).and.QDO(QTE)) M(id(QTT),id(QTE)) = Il(1,l)*Al(QTT,l)*Al(QTE,l)
    !if (QDO(QTT).and.QDO(QEE)) M(id(QTT),id(QEE)) = Il(2,l)*Al(QTT,l)*Al(QEE,l)
    !if (QDO(QTE).and.QDO(QEE)) M(id(QTE),id(QEE)) = Il(3,l)*Al(QTE,l)*Al(QEE,l)
    !if (QDO(QTB).and.QDO(QEB)) M(id(QTB),id(QEB)) = Il(4,l)*Al(QTB,l)*Al(QEB,l)
    !do X = 1, qn
    !  if (QDO(X)) M(id(X),id(X)) = Al(X,l)
    !  do Y = X + 1, qn
    !    if(QDO(X).and.QDO(Y)) M(id(Y),id(X)) = M(id(X),id(Y))
    !  end do
    !end do 
    !call inv_lapack(M)

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
    !do X = 1, qn
    !  if(QDO(X)) Nl(X,l) = sum(M(id(X),:))
    !end do

    deallocate(M)

  end do

end subroutine qmv


subroutine qall(QDO,lmax,rlmin,rlmax,fC,OC,Ag,Ac,Nlg,Nlc,gtype)
!*  Compute MV estimator normalization. Currently BB is ignored. 
!*
!*  Args:
!*    :QDO[6] (bool): Specifying which estimators to be combined for the minimum variance estimator, with size (6). The oder is TT, TE, EE, TB, EB and BB. 
!*    :lmax (int):    Maximum multipole of the output power spectra
!*    :rlmin/rlmax (int)   : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC/OC [l] (double): Theory/Observed CMB angular power spectra (TT, EE, BB, TE), with bounds (0:rlmax) 
!*
!*  Args(optional):
!*    :gtype (str): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [6,l] (double)  : Normalization of the TT, TE, EE, TB, EB, and MV estimators for lensing potential, with bounds (6,0:lmax)
!*    :Ac [6,l] (double)  : Same as Ag but for curl mode
!*    :Nlg [6,l] (double) : Weights for TT, TE, EE, TB, EB, and BB (=0) estimators for lensing potential, with bounds (6,0:lmax)
!*    :Nlc [6,l] (double) : Same as Nlg but for curl mode
!*
  implicit none
  !I/O
  logical, intent(in), dimension(6) :: QDO
  integer, intent(in) :: rlmin, rlmax, lmax
  double precision, intent(in), dimension(4,0:rlmax) :: fC, OC
  double precision, intent(out), dimension(6,0:lmax) :: Ag, Ac, Nlg, Nlc
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  character(1) :: gt
  integer :: TT = 1, EE = 2, BB = 3, TE = 4
  double precision, dimension(:,:), allocatable :: Ilg, Ilc

  gt = ''
  if (present(gtype)) gt = gtype

  !//// interface ////!
  Ag  = 0d0
  Ac  = 0d0
  Nlg = 0d0
  Nlc = 0d0
  if (QDO(1))  call qtt(lmax,rlmin,rlmax,fC(TT,:),OC(TT,:),Ag(1,:),Ac(1,:),gtype=gt)
  if (QDO(2))  call qte(lmax,rlmin,rlmax,fC(TE,:),OC(TT,:),OC(EE,:),Ag(2,:),Ac(2,:),gtype=gt)
  if (QDO(3))  call qee(lmax,rlmin,rlmax,fC(EE,:),OC(EE,:),Ag(3,:),Ac(3,:),gtype=gt)
  if (QDO(4))  call qtb(lmax,rlmin,rlmax,fC(TE,:),OC(TT,:),OC(BB,:),Ag(4,:),Ac(4,:),gtype=gt)
  if (QDO(5))  call qeb(lmax,rlmin,rlmax,fC(EE,:),OC(EE,:),OC(BB,:),Ag(5,:),Ac(5,:),gtype=gt)

  allocate(Ilg(4,0:lmax),Ilc(4,0:lmax))
  if (QDO(1).and.QDO(2))  call qttte(lmax,rlmin,rlmax,fC(TT,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(1,:),Ilc(1,:),gtype=gt)
  if (QDO(1).and.QDO(3))  call qttee(lmax,rlmin,rlmax,fC(TT,:),fC(EE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(2,:),Ilc(2,:),gtype=gt)
  if (QDO(2).and.QDO(3))  call qteee(lmax,rlmin,rlmax,fC(EE,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:),Ilg(3,:),Ilc(3,:),gtype=gt)
  if (QDO(4).and.QDO(5))  call qtbeb(lmax,rlmin,rlmax,fC(EE,:),fC(BB,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(BB,:),OC(TE,:),Ilg(4,:),Ilc(4,:),gtype=gt)
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
!*    :conv (double) : a parameter for convergence the iteration, default to 0.001
!*
!*  Returns:
!*    :Ag [l] (double) : CMB lensing potential normalization, with bounds (0:lmax)
!*    :Ac [l] (double) : Curl mode (pseudo lensing potential) normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, elmax, rlmin, rlmax, dlmin, dlmax
  double precision, intent(in), dimension(0:elmax) :: CE, OCE
  double precision, intent(in), dimension(0:rlmax) :: OCB
  double precision, intent(in), dimension(0:dlmax) :: Cpp
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  integer, intent(in), optional :: iter
  double precision, intent(in), optional :: conv
  !f2py integer :: iter = 1
  !f2py double precision :: conv = 1e-6
  !internal
  integer :: i, n, l, it
  double precision :: ratio, c
  double precision :: AgEB(0:dlmax), rCBB(0:rlmax)

  AgEB=0d0
  rCBB=0d0

  if (elmax<dlmax.or.elmax<rlmax) stop 'error (qeb_iter): does not support elmax<dlmax or elmax<rlmax'

  !initial values
  ratio = 1d0
  rCBB  = OCB
  it = 1
  c  = 0.001
  if (present(iter)) it = iter
  if (present(conv)) c  = conv

  do n = 1, it !loop for iteration 

    !lensing reconstruction with EB
    call qeb(dlmax,rlmin,rlmax,CE(0:rlmax),OCE(0:rlmax),rCBB,AgEB,Ac)

    !convergence check using gradient mode
    if (n>=2) then
      ratio = (sum(Ag)/sum(AgEB)-1d0)/dble(dlmax)
      write(*,*) n, ratio
    end if
    Ag = AgEB

    if (abs(ratio) < c) exit

    !delensing with EB-estimator
    call clbb_est((/rlmin,rlmax/),(/dlmin,dlmax/),CE(1:dlmax),Cpp(1:dlmax),OCE(1:dlmax)-CE(1:dlmax),AgEB,rCBB)
    rCBB = OCB - rCBB !delensed B-mode

    if(n==it) stop 'not converged'

  end do


end subroutine qeb_iter


subroutine ttt(lmax,rlmin,rlmax,fC,OCT,Ag,gtype)
!*  Cross normalization between lensing potential and amplitude modulation from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Args(optional):
!*    :gtype (str)       : Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
!*
!*  Returns:
!*    :Ag [l] (double)   : Cross normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag
  !optional
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: l, rL(2)
  double precision, dimension(lmax) :: lk2
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, G0

  write(*,*) 'norm tTT (lens)'
  rL = (/rlmin,rlmax/)

  lk2 = 1d0
  if (present(gtype).and.gtype=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)
    end do
  end if

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_lens.ttt): observed cltt is zero'
  end do

  W1 = 1d0 / OCT(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  S0 = 0d0
  call kernels_lenstau(rL,W1,W2,S0,'S0')

  W2 = W1 * fC(rlmin:rlmax)
  G0 = 0d0
  call kernels_lenstau(rL,W2,W2,G0,'G0')

  Ag = 0d0
  do l = 1, lmax
    Ag(l) = (S0(1,l)+G0(1,l))/lk2(l)
  end do

end subroutine ttt


end module norm_lens


