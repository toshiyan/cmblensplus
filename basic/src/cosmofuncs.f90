!////////////////////////////////////////////////////!
! * Functions for Cosmology
!////////////////////////////////////////////////////!

module cosmofuncs
  use constants, only: pi
  use cosmofunc, only: C_z, H_z, dL_dz, dH_dz, D_z, g_factor, g_rate, cosmoparams, zoftau
  use utilsgw,   only: dotn, dlndotn
  implicit none

  private C_z, H_z, dL_dz, dH_dz, D_z, g_factor

contains


subroutine hubble(z,H0,Om,Ov,w0,wa,zn,Hz)
!*  Compute the expansion rate in unit of 1/Mpc, H/c. 
!*
!*  Args:
!*    :z[zn] (double)   : Redshifts at which H is computed
!*
!*  Args(optional):
!*    :H0 (double)      : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)      : The current value of Omega_matter, default to 0.3
!*    :Ov (double)      : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double)  : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*
!*  Returns:
!*    :Hz[zn] (double)   : The expansion rate, H(z)/c (H is divided by c). 
!*

  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: Hz
  integer :: j
  type(cosmoparams) :: cp
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    Hz(j) = H_z(z(j),cp)
  end do

end subroutine hubble


subroutine dhubble_dz(z,H0,Om,Ov,w0,wa,zn,dHdz)
  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: dHdz
  integer :: j
  type(cosmoparams) :: cp
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    dHdz(j) = dH_dz(z(j),cp)
  end do

end subroutine dhubble_dz


subroutine dist2z(rz,H0,Om,Ov,w0,wa,zn,z)
  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: rz
  double precision, intent(out), dimension(0:zn-1) :: z
  integer :: j
  type(cosmoparams) :: cp
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(rz)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    z(j) = zoftau(-rz(j),0d0,cp)
  end do

end subroutine dist2z


subroutine dist_comoving(z,H0,Om,Ov,w0,wa,zn,rz)
  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: rz
  integer :: j
  type(cosmoparams) :: cp
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    rz(j) = C_z(z(j),cp)
  end do

end subroutine dist_comoving


subroutine dist_luminosity(z,H0,Om,Ov,w0,wa,zn,DLz)
  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: DLz
  integer :: j
  type(cosmoparams) :: cp
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    DLz(j) = (1+z(j))*C_z(z(j),cp)
  end do

end subroutine dist_luminosity


subroutine growth_factor(z,H0,Om,Ov,w0,wa,zn,Dz,normed)
  implicit none
  logical, intent(in) :: normed
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: Dz
  integer :: j
  double precision :: wz
  type(cosmoparams) :: cp
  !opt4py :: normed = False
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    if (normed) then
      Dz(j) = D_z(z(j),cp)
    else
      wz = cp%w0 + (1d0-1d0/(1d0+z(j)))*cp%wa
      Dz(j) = g_factor(z(j),cp%Om,cp%Ov,wz)
    end if
  end do

end subroutine growth_factor


subroutine growth_rate(z,H0,Om,Ov,w0,wa,zn,fz)
  implicit none
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: fz
  integer :: j
  type(cosmoparams) :: cp
  !opt4py :: H0 = 70.
  !opt4py :: Ov = 0.7
  !opt4py :: Om = 0.3
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    fz(j) = g_rate(z(j),cp)
  end do

end subroutine growth_rate



subroutine nz_gw(z,Cz,Hz,ntype,dotn0,Tobs,nz)
!*  Distribution function of NS binary sources per redshift (dN/dz)
!
  implicit none
! [inputs]  
!   ntype  --- type of dotn functional from
!   z      --- redshift
!   Cz     --- comoving distance
!   Hz     --- expansion rate
  character(4), intent(in) :: ntype
  double precision, intent(in) :: z, Cz, Hz
! (optional)
!   dotn0   --- current merger-rate
!   Tobs    --- total observation time
  double precision, intent(in) :: dotn0, Tobs
! [output]
!   f       --- distrubution function at z
  double precision, intent(out) :: nz
! [internal]
  double precision :: Lz, sz
  !opt4py :: dotn0 = 1e-6
  !opt4py :: Tobs = 3.

  Lz = Cz*(1d0+z)
  sz = dotn(z,ntype)
  nz = Tobs*dotn0*4d0*pi*sz*Lz**2/(Hz*(1d0+z)**3)

end subroutine nz_gw


subroutine drate_dz(z,zn,ntype,dndz)
  implicit none
  character(4), intent(in) :: ntype
  integer, intent(in) :: zn
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: dndz
  integer :: j
  !opt4py :: ntype = 'CH06'
  !opt4py :: zn = 0
  !add2py :: zn = len(z)

  do j = 0, zn-1
    dndz(j) = dlndotn(z(j),ntype)
  end do

end subroutine drate_dz


end module cosmofuncs


