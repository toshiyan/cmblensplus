!/////////////////////////////////////////////////////////////////!
! * Calculation in flatsky 
! - Not recommended for actual calculations, just for cross check
!/////////////////////////////////////////////////////////////////!

module flat
  use norm_flat
  use bb_flat
  implicit none

contains


subroutine alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,Ag,Ac,gln,gle,lxcut)
!*  Compute flat-sky quadratic estimator normalization
!*  CAUTION: This code interpolates the input Cl at the non-integer multipole by simply Cl(int(ell)) which leads to a small discrepancy in the normalization computed from the FFT-based method (which uses linear interpolation) and from this code. It is desireble to use the FFT-based normalization if you want to normalize the simulation results.
!* 
!*  Args:
!*    :qest (str)  : estimator combination (TT, TE, TB, EE, EB, or BB)
!*    :qtype (str) : estimator type (lensing, patchytau)
!*    :lmax (double) : output maximum multipole
!*    :rlmax/rlmin (double) : input CMB multipole range for reconstruction
!*    :fC[rlmax] (double) : power spectrum in the numerator
!*    :W1/W2[rlmax] : inverse of the observed power spectrum
!*
!*  Args(optional):
!*    :gln (int) : number of the GL integration points
!*    :lxcut (int) : multipole cut in x-direction, |l_x| < lx
!*    :gle (double) : convergence parameter for the GL integration
!*
!*  Returns:
!*    :Ag/Ac[l] (double) : normalization for even and odd estimator pairs
!*
  !f2py intent(in) qest, qtype, lmax, rlmin, rlmax, fC, W1, W2, gln, lxcut, gle
  !f2py intent(out) Ag, Ac
  !f2py depend(rlmax) fC, W1, W2
  !f2py depend(lmax) Ag, Ac
  implicit none
  !I/O
  character(*), intent(in) :: qest, qtype
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, W1, W2
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  integer, intent(in), optional :: gln, lxcut
  double precision, intent(in), optional :: gle
  !f2py integer :: gln = 100
  !f2py integer :: lxcut = 0
  !f2py double precision :: gle = 1e-14

  Ag(0) = 0d0
  Ac(0) = 0d0
  call alxy_flat_integ(qest,qtype,(/1,lmax/),(/rlmin,rlmax/),Ag(1:),Ac(1:),fC(1:),W1(1:),W2(1:),gln=gln,gle=gle,lxcut=lxcut)

end subroutine alxy

subroutine alxy_asym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,Ag,Ac,gln,gle,lxcut)
  !f2py intent(in) qest, qtype, lmax, rlmin, rlmax, fC, AA, BB, AB, gln, lxcut, gle
  !f2py intent(out) Ag, Ac
  !f2py depend(rlmax) fC, AA, BB, AB
  !f2py depend(lmax) Ag, Ac
  implicit none
  !I/O
  character(*), intent(in) :: qest, qtype
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, AA, BB, AB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  integer, intent(in), optional :: gln, lxcut
  double precision, intent(in), optional :: gle
  !f2py integer :: gln = 100
  !f2py integer :: lxcut = 0
  !f2py double precision :: gle = 1e-14

  Ag(0) = 0d0
  Ac(0) = 0d0
  call alxy_flat_integ(qest,qtype,(/1,lmax/),(/rlmin,rlmax/),Ag(1:),Ac(1:),fC(1:),AA=AA(1:),BB=BB(1:),AB=AB(1:),gln=gln,gle=gle,lxcut=lxcut)

end subroutine alxy_asym

subroutine bbxy(lmax,rlmin,rlmax,XX,YY,BB,weight,gln,gle)
  !f2py intent(in) lmax, rlmin, rlmax, XX, YY, weight, gln, gle
  !f2py intent(out) BB
  !f2py depend(rlmax) XX, YY
  !f2py depend(lmax) BB
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: XX, YY
  double precision, intent(out), dimension(2,0:lmax) :: BB
  !optional
  character(*), intent(in), optional :: weight
  integer, intent(in), optional :: gln
  double precision, intent(in), optional :: gle
  !f2py character(*) :: weight = 'lensing'
  !f2py integer :: gln = 100
  !f2py double precision :: gle = 1e-14

  call bbxy_flat((/rlmin,rlmax/),(/1,lmax/),BB,XX,YY,weight,gln,gle)

end subroutine bbxy



end module flat



