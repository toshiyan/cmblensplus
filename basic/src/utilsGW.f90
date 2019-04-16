!////////////////////////////////////////////////////!
! * GW detector forecast
!////////////////////////////////////////////////////!

module utilsGW
  use constants, only: pi
  use cosmofunc, only: C_z, H_z, dL_dz, dH_dz, cosmoparams
  use general,   only: linspace
  implicit none

  private C_z, H_z, dL_dz, dH_dz, linspace

contains


function mean_Lz(zb,cp,cp0,ntype)  result(f)
! * Mean luminosity-distance times total number (d_L x N_i)
!
  implicit none
  type(cosmoparams), intent(in) :: cp, cp0
  integer, intent(in) :: ntype
  double precision, intent(in) :: zb(1:2)
  integer :: j
  double precision :: f, z(1:10000), nd, dz, CzHz(1:2), CzHzf(1:2)

  f  = 0d0
  z  = linspace(zb,10000)
  dz = z(2)-z(1)
  do j = 1, 10000
    CzHzf = [C_z(z(j),cp0), H_z(z(j),cp0)]
    nd    = nz_gw(z(j),CzHzf,ntype)/(dL_dz(z(j),CzHzf)) ! fixed cosmology
    ! changing variable from d_L to z
    CzHz  = [C_z(z(j),cp), H_z(z(j),cp)]
    f     = f + dz*dL_dz(z(j),CzHz) * CzHz(1)*(1d0+z(j)) * nd 
  end do

end function mean_Lz


function mean_Nz(zb,cp,cp0,ntype)  result(f)
! * Total number (N_i)
!
  implicit none
  type(cosmoparams), intent(in) :: cp, cp0
  integer, intent(in) :: ntype
  double precision, intent(in) :: zb(1:2)
  integer :: j
  double precision :: f, z(1:10000), nd, dz, CzHz(1:2), CzHzf(1:2)

  f  = 0d0
  z  = linspace(zb,10000)
  dz = z(2)-z(1)
  do j = 1, 10000
    CzHzf = [C_z(z(j),cp0), H_z(z(j),cp0)]
    nd    = nz_gw(z(j),CzHzf,ntype)/(dL_dz(z(j),CzHzf)) ! fixed cosmology
    ! changing variable from d_L to z
    CzHz  = [C_z(z(j),cp), H_z(z(j),cp)]
    f     = f + dz*dL_dz(z(j),CzHz) * nd 
  end do

end function mean_Nz


function nz_gw(z,CzHz,ntype,dotn0,Tobs)  result(f)
! * Distribution function of NS binary sources per redshift (dN/dz)
!
  implicit none
! [inputs]  
!   ntype   --- type of dotn functional from
!   z       --- redshift
!   CzHz(2) --- array containing [comoving distance, expansion rate]
  integer, intent(in) :: ntype
  double precision, intent(in) :: z, CzHz(1:2)
!
! (optional)
!   dotn0   --- current NS-NS merger-rate
!   Tobs    --- total observation time
  double precision, intent(in), optional :: dotn0, Tobs
!
! [output]
!   f       --- distrubution function at z
  double precision :: f
!
! [internal]
  double precision :: Lz, Hz

  !* set luminosity distance and expansion rate
  Lz = CzHz(1)*(1d0+z)
  Hz = CzHz(2)

  !* distribution function
  f = 4d0*pi*dotn(z,ntype)*Lz**2/(Hz*(1d0+z)**3)

  !* corrections
  if (present(dotn0)) f = f*dotn0
  if (present(Tobs))  f = f*Tobs

end function nz_gw


function Fz_gw(z,CzHz,CzHzf,dHdzf)  result(f)
! * correction factor to the distribution function (Fz_gw = K(z)/n_d)
!
  implicit none
! [inputs]  
!   z        --- redshift
!   CzHz(2)  --- array containing [comoving distance, expansion rate]
!   CzHzf(2) --- array containing [comoving distance, expansion rate]
  double precision, intent(in) :: z, CzHz(1:2), CzHzf(1:2), dHdzf
!
! [output]
!   f           --- distrubution function at z
  double precision :: f
!
! [internal]
  double precision :: a, Lz, Lzf, Hzf, dLdz, dLdzf, n_d

  !* set luminosity distance and expansion rate
  Lzf   = CzHzf(1)*(1d0+z)
  Hzf   = CzHzf(2)
  dLdzf = dL_dz(z,CzHzf)

  !* calculating dln(n_d)/dz (fixed cosmology)
  f = 0d0
  if (z>0d0) f = 2d0*(dLdzf/Lzf - 1d0/(1d0+z)) - (2d0*(1d0+z)+dLdzf*Hzf+Lzf*dHdzf)/((1d0+z)**2+Lzf*Hzf)
  if (dotn(z)>0d0)  f = f + ddotn(z)/dotn(z)

  !* to dln(n_d)/dD = dln(n_d)/dz * dz/dD (fixed cosmology)
  f = f / dLdzf

  !* make correction factor ( dD/dz * D * [ 1 + D * dln(n_d)/dD ] )
  Lz   = CzHz(1)*(1d0+z)
  dLdz = dL_dz(z,CzHz)
  f    = dLdz * Lz * ( 1d0 + Lz * f)

end function Fz_gw


function dotn(z,ntype)  result(f)
! * redshift dependence of NS-NS merger rate
!
  implicit none
  integer, intent(in), optional :: ntype
  double precision, intent(in) :: z
  integer :: n
  double precision :: f

  n = 0
  if (present(ntype)) n=ntype

  select case (n)
  case(0) !Cutler & Harmus (2006) fitting formula
    if (z<=1d0) then
      f = 1d0+2d0*z
    else if (z>1d0.and.z<=5d0) then 
      f = 0.75d0*(5d0-z)
    else 
      f = 0d0
    end if
  case(1) 
    f = 1d0
  end select

end function dotn


function ddotn(zs,zb,n0)  result(f)
  implicit none
  double precision, intent(in) :: zs
  double precision, intent(in), optional :: zb(1:2), n0
  double precision :: f

  if (zs<=1d0) then
    f = 2d0
  else if (zs>1d0.and.zs<=5d0) then 
    f = -0.75d0
  else 
    f = 0d0
  end if

  if (present(zb)) then
    if (zs<zb(1).or.zb(2)<zs) f=0d0
  end if
  if (present(n0))  f = f*n0

end function ddotn


end module utilsGW

