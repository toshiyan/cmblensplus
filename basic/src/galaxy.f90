!/////////////////////////////////////////////////////////////////////!
! Some useful functions for computing galaxy-survey related quantities
!/////////////////////////////////////////////////////////////////////!

module galaxy
  use utilsgal, only: zbin_SF, ngal_SF, nz_SF_scal
  use funcs,    only: lnGamma
  implicit none

contains


subroutine dndz_sf(zn,z,a,b,zm,dndz)
!* Galaxy z distribution
!*
!*  Args:
!*    :z[zn] (double) : redshifts at which dNdz is returned
!*    :a, b (double)  : shape parameters of Schechter-like galaxy distribution
!*    :zm (double)    : mean redshift
!*
!*  Returns:
!*    :dndz[zn] (double) : galaxy z distribution
!*
  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: a, b, zm
  double precision, intent(in), dimension(1:zn) :: z
  double precision, intent(out), dimension(1:zn) :: dndz
  integer :: i
  !opt4py :: zn = 0
  !add2py :: if zn==0: zn=len(z)

  do i = 1, zn
    dndz(i) = nz_SF_scal(z(i),a,b,zm)
  end do

end subroutine dndz_sf


subroutine zbin(zn,a,b,zm,zb,verbose)
!* Computing z-interval of z-bin so that number of galaxies at each z-bin is equal
!*
!*  Args:
!*    :zn (int)       : number of z-bins
!*    :a, b (double)  : shape parameters of Schechter-like galaxy distribution
!*    :zm (double)    : mean redshift
!*
!*  Args(optional):
!*    :verbose (bool) : output messages (default to True)
!*
!*  Returns:
!*    :zb[zn+1] (double) : z-intervals
!*
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: zn
  double precision, intent(in) :: a, b, zm
  double precision, intent(out), dimension(1:zn+1) :: zb
  !opt4py :: verbose = False

  call zbin_SF(a,b,zm,zb)

  if (verbose)  write(*,*) 'zbin =', zb

end subroutine zbin


subroutine frac(zn,zb,a,b,zm,nfrac,zbias,sigma,verbose)
!* Computing z-interval of z-bin so that number of galaxies at each z-bin is equal
!*
!*  Args:
!*    :zn (int)       : number of z-bins
!*    :a, b (double)  : shape parameters of Schechter-like galaxy distribution
!*    :zm (double)    : mean redshift
!*    :zb[zn+1] (double) : z-intervals
!*
!*  Args(optional):
!*    :verbose (bool) : output messages (default to True)
!*    :zbias (double) : constant bias to true photo-z
!*    :sigma (double) : uncertaines of photo-z
!*
!*  Returns:
!*    :nfrac[zn] (double) : fraction of galaxy number at each bin, defined by int_zi^zi+1 dz N(z)/int dz N(z) 
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: zn
  double precision, intent(in) :: a, b, zm, sigma, zbias
  double precision, intent(in), dimension(1:zn+1) :: zb
  double precision, intent(out), dimension(1:zn) :: nfrac
  !opt4py :: verbose = False
  !opt4py :: zbias = 0.0
  !opt4py :: sigma = 0.0

  call ngal_SF(nfrac,zb,a,b,zm,sigma,zbias)
  
  if (verbose)  write(*,*) 'fractional Ngal = ', nfrac

end subroutine frac


end module galaxy


