!////////////////////////////////////////////////////!
! * Angular power spectrum calculation
!////////////////////////////////////////////////////!

module aps
  use constants, only: pi
  use pstool, only: binned_ells
  implicit none

  private pi
  private binned_ells

contains


subroutine binning(bn,eL,bp,bc,spc)
!* Return multipole-bin edges and centers
!* 
!* Args:
!*   :bn (int)    : number of bins
!*   :eL[2] (int) : bin edges
!*
!* Args(optional):
!*   :spc (str)   : bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing
!*
!* Returns:
!*   :bp (double) : bin edges, with bounds (0:bn)
!*   :bc (double) : bin centers, with bounds (bn)
!*
  implicit none
  !I/O
  character(4), intent(in) :: spc
  integer, intent(in) :: bn
  integer, intent(in), dimension(2) :: eL
  double precision, intent(out), dimension(0:bn) :: bp
  double precision, intent(out), dimension(bn) :: bc
  !optional
  !opt4py :: spc = ''

  call binned_ells(eL,bp,bc,spc=spc)

end subroutine binning


subroutine map_vars(lmax,cl,sigma)
!*  Variance of a map and its derivative
!* 
!*  Args:
!*    :lmax (int)    : Maximum multipole of the input cl
!*    :cl[l] (double): Angular power spectrum, with bounds (0:lmax)
!*
!*  Returns:
!*    :sigma[2] (double): Variance of the map and of map derivative, with bounds (2)
!*  
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in), dimension(0:lmax) :: cl
  double precision, intent(out), dimension(0:1) :: sigma
  !internal
  integer :: l
  double precision :: al
  
  sigma = 0d0
  do l = 1, lmax
    al = dble(l)
    sigma(0) = sigma(0) + (2d0*al+1d0)*cl(l)
    sigma(1) = sigma(1) + (2d0*al+1d0)*(al**2+al)*cl(l)
  end do 
  sigma = sigma/(4d0*pi)

end subroutine map_vars


end module aps

