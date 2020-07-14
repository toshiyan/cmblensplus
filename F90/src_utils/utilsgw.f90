!////////////////////////////////////////////////////!
! * Functions for Cosmology
!////////////////////////////////////////////////////!

module utilsgw
  use constants, only: pi
  use cosmofunc, only: C_z, H_z, dL_dz, dH_dz, cosmoparams
  implicit none

  private C_z, H_z, dL_dz, dH_dz

contains


function dotn(z,ntype) result(f)
!*  Redshift dependence of NS-NS merger rate
  implicit none
  character(4), intent(in) :: ntype
  double precision, intent(in) :: z
  double precision :: f

  select case (ntype)
  case('CH06') !Cutler & Harmus (2006) fitting formula
    if (z<=1d0) then
      f = 1d0+2d0*z
    else if (z>1d0.and.z<=5d0) then 
      f = 0.75d0*(5d0-z)
    else 
      f = 0d0
    end if
  case('none') 
    f = 1d0
  case default
    write(*,*) 'dotn: unknown z-distribution type of gw sources', ntype
    stop 
  end select

end function dotn


function ddotn(zs,ntype) result(f)
  implicit none
  character(4), intent(in) :: ntype
  double precision, intent(in) :: zs
  double precision :: f

  select case (ntype)
  case('CH06') !Cutler & Harmus (2006) fitting formula
    if (zs<=1d0) then
      f = 2d0
    else if (zs>1d0.and.zs<=5d0) then 
      f = -0.75d0
    else 
      f = 0d0
    end if
  case('none')
    f = 0d0
  case default
    write(*,*) 'dotn: unknown z-distribution type of gw sources', ntype
    stop
  end select

end function ddotn


function dlndotn(zs,ntype) result(f)
  implicit none
  character(4), intent(in) :: ntype
  double precision, intent(in) :: zs
  double precision :: f, dsz, sz

  sz  = dotn(zs,ntype)
  dsz = ddotn(zs,ntype)
  if (sz==0d0) then
    f = 0d0
  else
    f = dsz/sz
  end if

end function dlndotn


end module utilsgw


