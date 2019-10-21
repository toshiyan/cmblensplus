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


