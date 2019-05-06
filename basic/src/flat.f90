!/////////////////////////////////////////////////////////////////!
! * Calculation in flatsky 
! - Not recommended for actual calculations, just for cross check
!/////////////////////////////////////////////////////////////////!

module flat
  use norm_flat
  use bb_flat
  implicit none

contains


subroutine alxy(est,lmax,rlmin,rlmax,fC,W1,W2,Ag,Ac,qe,gln,gle,lxcut)
  implicit none
  !I/O
  character(*), intent(in) :: est
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, W1, W2
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !optional
  character(8), intent(in), optional :: qe
  integer, intent(in), optional :: gln, lxcut
  double precision, intent(in), optional :: gle
  !f2py character(8) :: qe = 'lensing'
  !f2py integer :: gln = 100
  !f2py integer :: lxcut = 0
  !f2py double precision :: gle = 1e-14

  Ag(0) = 0d0
  Ac(0) = 0d0
  call alxy_flat(est,(/1,lmax/),(/rlmin,rlmax/),Ag(1:),Ac(1:),fC(1:),W1(1:),W2(1:),qe,gln,gle,lxcut)

end subroutine alxy


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


