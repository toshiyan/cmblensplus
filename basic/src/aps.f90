!////////////////////////////////////////////////////!
! * Angular power spectrum calculation
!////////////////////////////////////////////////////!

module aps
  use constants, only: pi
  use pstool, only: binned_ells, readcl_camb, power_binning
  implicit none

  private pi
  private binned_ells, readcl_camb, power_binning

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
  integer, intent(in) :: bn
  integer, intent(in), dimension(2) :: eL
  double precision, intent(out), dimension(bn) :: bc
  double precision, intent(out), dimension(0:bn) :: bp
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
!*  Return CMB cls from CAMB output files
!*
!*  Args:
!*    :f (str): Filename 
!*    :lmin (int)      : Minimum multipole of the output Cl
!*    :lmax (int)      : Maximum multipole of the output Cl
!*    :numcls (int)    : Number of cls to be read
!*
!*  Args(Optional):
!*    :bb (bool)   : Including BB in the file or not. The data should be TT, EE, TE, dd, Td (,Ed) if bb = False (default), and TT, EE, BB, TE if bb = True. 
!*    :raw (bool)  : The cls in the file are multiplied by l(l+1)/2pi if raw = False (default)
!*
!*  Returns:
!*    :cl[numcls,l] (double): Angular power spectra, with bountd (numcls,0:lmax)
!*
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


subroutine cl2bcl(bn,lmax,cl,cb,spc)
!*  From unbinned to binned angular power spectrum
!*
!*  Args:
!*    :bn (int)         : number of multipole bins
!*    :lmax (int)       : maximum multipole of the input angular power spectrum
!*    :cl[l] (double)   : angular power spectrum, with bounds (0:lmax)
!*
!* Args(optional):
!*    :spc (str) : bin spacing, '' = linear (default), 'log' = log spacing, 'log10' = log10 spacing, 'p2' = power of 2 spacing, 'p3' = power of 3 spacing
!*
!*  Returns:
!*    :cb[bin] (double) : auto or cross angular power spectrum with multipole binning, with bounds (0:bn-1)
!*
  implicit none
  integer, intent(in) :: bn, lmax
  double precision, intent(in), dimension(0:lmax) :: cl
  double precision, intent(out), dimension(bn) :: cb
  !optional
  character(4), intent(in), optional :: spc
  !f2py character(4) :: spc = ''
  !internal
  integer :: eL(2)
  double precision, allocatable :: bp(:)
  character(4) :: sp

  sp = ''
  if (present(spc)) sp = spc

  eL = (/1,lmax/)

  allocate(bp(bn+1))
  call binned_ells(eL,bp,spc=sp)
  call power_binning(bp,eL,cl(1:lmax),Cb)
  deallocate(bp)

end subroutine cl2bcl


end module aps

