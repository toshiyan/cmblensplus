!////////////////////////////////////////////////////!
! * Galaxy survey forecast
!////////////////////////////////////////////////////!

module utilsgal
  use funcs,   only: erfc, lnGamma
  use general, only: linspace
  implicit none

  private erfc, lnGamma
  private linspace

contains


subroutine z02zm(a,b,z0,zm)
!*  Transform z0 to zm for Shechter-like galaxy distribution
!*  
!*  Args:
!*    :a, b (double) : shape parameters of Schechter-like galaxy distribution
!*    :z0 (double)   : a characteristic redshift of Schechter-like galaxy distribution
!*
!*  Returns:
!*    :zm (double) : mean redshift
!*
  implicit none
  double precision, intent(in) :: a, b, z0
  double precision, intent(out) :: zm

  zm = z0*dexp(lnGamma((a+2d0)/b)-lnGamma((a+1d0)/b))

end subroutine z02zm


subroutine zm2z0(a,b,zm,z0)
!*  Transform zm to z0 for Shechter-like galaxy distribution
!*  
!*  Args:
!*    :a, b (double) : shape parameters of Schechter-like galaxy distribution
!*    :zm (double)   : a characteristic redshift of Schechter-like galaxy distribution
!*
!*  Returns:
!*    :z0 (double) : mean redshift
!*
  implicit none
  double precision, intent(in) :: a, b, zm
  double precision, intent(out) :: z0

  z0 = zm*dexp(lnGamma((a+1d0)/b)-lnGamma((a+2d0)/b))

end subroutine zm2z0


function nz_SF_scal(z,a,b,zm)  result(f)
  implicit none
  double precision, intent(in) :: z, a, b, zm
  double precision :: N, z0, f

  call zm2z0(a,b,zm,z0)
  N = b / (z0*dexp(lnGamma((a+1d0)/b)))
  
  f = N*(z/z0)**a*dexp(-(z/z0)**b)

end function nz_SF_scal


function nz_SF_arr(z,a,b,zm) result(f)
  implicit none
  double precision, intent(in) :: z(:), a, b, zm
  integer :: i
  double precision :: N, z0, f(size(z))

  !call zm2z0(a,b,zm,z0)
  !N = b / (z0*dexp(lnGamma((a+1d0)/b)))

  do i = 1, size(z)
    f(i) = nz_SF_scal(z(i),a,b,zm)
    !f(i) = N*(z(i)/z0)**a*dexp(-(z(i)/z0)**b)
  end do

end function nz_SF_arr


function nz_21cm(z,nzp) result(f)
! - Oct 15, 2015
  implicit none
  double precision, intent(in) :: z, nzp(1:3)
  double precision :: f

  f = z**(nzp(1))*dexp(-nzp(2)*z**nzp(3))

end function nz_21cm


function nz_21cm_histgram(z,fluxcut) result(f)
! * older version
  implicit none
  character(*), intent(in) :: fluxcut
  double precision, intent(in) :: z
  double precision :: N(12), f

  select case(fluxcut)
  case("100nJy")
    N = [10815d0,25483d0,27587d0,23100d0,21093d0,16764d0,14190d0,12646d0,9927d0,8462d0,7405d0,6396d0]
  case("1uJy")
    N = [5445d0,10987d0,12171d0,10327d0,7655d0,5434d0,3870d0,2928d0,2216d0,1675d0,1344d0,1076d0]
  case("5uJy")
    N = [2642d0,4499d0,4405d0,3381d0,2155d0,1403d0,932d0,648d0,448d0,317d0,231d0,172d0]
  case("10uJy")
    N = [1701d0,2731d0,2548d0,1875d0,1126d0,690d0,431d0,282d0,186d0,125d0,88d0,65d0]
  end select

  f = N(int(2*z)+1)/sum(N)

end function nz_21cm_histgram


function nz_delta_scal(z,zmean,zwidth) result(f)
  implicit none
  double precision, intent(in) :: z, zmean, zwidth
  double precision :: f

  f = 0d0
  if(abs(z-zmean)<zwidth)  f = 1d0

end function nz_delta_scal


function pz_SF_scal(z,zi,sigma,zbias)  result(f)
  implicit none
  double precision, intent(in) :: z, zi(1:2), sigma, zbias
  double precision :: f, s, c

  f = 0d0
  c = zbias*(1d0+z)
  if(sigma==0d0) then
    if(z-c>=zi(1).and.z-c<zi(2)) f = 1d0
  else
    s = sigma*(1d0+z)*dsqrt(2d0)
    f = (erfc((zi(1)-z+c)/s)-erfc((zi(2)-z+c)/s))*0.5d0
  end if

end function pz_SF_scal


function pz_SF_arr(z,zi,sigma,zbias)  result(f)
  implicit none
  !I/O
  double precision, intent(in) :: z(:), zi(1:2), sigma, zbias
  !internal
  integer :: i
  double precision :: f(size(z)), s, c

  do i = 1, size(z)
    f(i) = pz_SF_scal(z(i),zi,sigma,zbias)
  end do

end function pz_SF_arr


subroutine zbin_SF(a,b,zm,zb)
!* divide total galaxy distribution (ns) into z-bins
  implicit none
  double precision, intent(in) :: a, b, zm
  double precision, intent(out) :: zb(:)
  integer :: i, j, n, nz
  integer, parameter :: jn = 5000
  double precision :: z(jn), g(jn), zz(size(zb)-1), gmax, zmax

  zmax = 20d0
  z = linspace(0d0,zmax,jn)
  nz = size(zb)-1

  gmax = sum(nz_SF_arr(z,a,b,zm))*(z(2)-z(1))
  do j = 1, jn
    z = linspace(0d0,j*zmax/dble(jn),jn)
    g(j) = sum(nz_SF_arr(z,a,b,zm))*(z(2)-z(1))
    do n = 1, nz
      if(g(j) < n*gmax/dble(nz)) zz(n) = j*zmax/dble(jn)
    end do
  end do

  zb(1) = 0d0
  zb(2:nz+1) = zz(1:nz)

end subroutine zbin_SF


function integ_SF(zb,a,b,zm,sigma,zbias)  result(f)
  implicit none
  !I/O
  double precision, intent(in) :: zb(1:2), a, b, zm, sigma, zbias
  !internal
  integer :: i
  double precision :: f, z(1:100000)

  z = linspace([0d0,10d0],100000)
  f = 0d0
  do i = 1, 100000
    f = f + nz_SF_scal(z(i),a,b,zm) * pz_SF_scal(z(i),zb,sigma,zbias)
  end do
  f = f * (z(2)-z(1))  !multiply dz

end function integ_SF


subroutine ngal_SF(frac,zb,a,b,zm,sigma,zbias)
  implicit none
  !I/O
  double precision, intent(in) :: zb(:), a, b, zm, sigma, zbias
  double precision, intent(out) :: frac(:)
  !internal
  integer :: n, i, num

  num  = size(frac)
  frac = 0d0
  if (size(zb)/=num+1) stop 'error (ngal_SF): size of zb is strange'

  do n = 1, num
    frac(n) = integ_SF(zb(n:n+1),a,b,zm,sigma,zbias)
  end do
  frac = frac/sum(frac)

end subroutine ngal_SF


function gbias(z,btype)  result(f)
  implicit none
  character(*), intent(in) :: btype
  double precision, intent(in) :: z
  double precision :: f

  select case(btype)
  case('sqrtz')
    f = dsqrt(1d0+z)
  case default
    f = 1d0
  end select

end function gbias


end module utilsgal


