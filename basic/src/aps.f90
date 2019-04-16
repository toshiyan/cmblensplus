!////////////////////////////////////////////////////!
! * Angular power spectrum calculation
!////////////////////////////////////////////////////!

module aps
  use pstool, only: binned_ells, readcl_camb
  implicit none

  private binned_ells, readcl_camb

contains


subroutine binning(bn,eL,bp,bc,spc)
! * return binned multipole edges and centers
  implicit none
  !I/O
  integer, intent(in) :: bn
  integer, intent(in), dimension(2) :: eL
  double precision, intent(out), dimension(bn) :: bc
  double precision, intent(out), dimension(bn+1) :: bp
  !optional
  character(4), intent(in), optional :: spc
  !f2py character(4) :: spc = ''
  !internal
  character(4) :: sp

  sp = ''
  if (present(spc)) sp = spc
  call binned_ells(eL,bp,bc,sp)

end subroutine binning


subroutine read_cambcls(f,lmin,lmax,numcls,cl,bb,raw)
! Return CMB cls from CAMB output files
!   if bb = False: TT, EE, TE, dd, Td(, Ed)
!   if bb = True:  TT, EE, BB, TE
! raw = True if cl is not multiplied by l(l+1)/2pi
  implicit none
  character(*), intent(in) :: f
  integer, intent(in) :: lmin, lmax, numcls
  double precision, intent(out), dimension(numcls,0:lmax) :: cl
  !optional
  logical, intent(in), optional :: bb, raw
  !f2py logical :: bb = 0
  !f2py logical :: raw = 0

  cl(:,0) = 0d0
  call readcl_camb(cl(:,1:lmax),f,(/lmin,lmax/),bb,raw)


end subroutine read_cambcls


end module aps

