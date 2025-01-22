!////////////////////////////////////////////////////!
! * Functions for Cosmology
!////////////////////////////////////////////////////!

module cosmofuncs
  use constants, only: pi, const_c
  use cosmofunc, only: C_z, H_z, dL_dz, dH_dz, D_z, g_factor, g_rate, cosmoparams, zoftau
  use utilsgw,   only: dotn, dlndotn
  implicit none

  private pi, const_c
  private C_z, H_z, dL_dz, dH_dz, D_z, g_factor, g_rate, cosmoparams, zoftau
  private dotn, dlndotn

contains


subroutine hubble(z,H0,Om,Ov,w0,wa,zn,Hz,divc)
!*  Compute the expansion rate in unit of 1/Mpc, H/c, or in unit of km/s/Mpc, H.  
!*
!*  Args:
!*    :z[zn] (double)  : Redshifts at which H is computed
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*    :divc (bool)     : Divide H by c or not, default to False. 
!*
!*  Returns:
!*    :Hz[zn] (double) : The expansion rate, H(z), divided by c or not. 
!*
  implicit none
  !f2py intent(in) zn, divc, H0, Om, Ov, w0, wa, z
  !f2py intent(out) Hz
  !f2py depend(zn) z, Hz
  integer, intent(in) :: zn
  logical, intent(in) :: divc
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: Hz
  integer :: j
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0., divc = False
  !add2py :: zn = len(z)

  cp%H0 = H0
  cp%Om = Om
  cp%Ov = Ov
  cp%w0 = w0
  cp%wa = wa

  do j = 0, zn-1
    Hz(j) = H_z(z(j),cp)
  end do

  if(.not.divc) Hz = Hz*const_c

end subroutine hubble

subroutine dhubble_dz(z,H0,Om,Ov,w0,wa,zn,dHdz)
!*  Compute dH(z)/dz.  
!*
!*  Args:
!*    :z[zn] (double)  : Redshifts at which dH/dz is computed
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*
!*  Returns:
!*    :dHdz[zn] (double) : The derivative of the expansion rate, dH(z)/dz. 
!*
  implicit none
  !f2py intent(in) zn, H0, Om, Ov, w0, wa, z
  !f2py intent(out) dHdz
  !f2py depend(zn) z, dHdz
  !I/O
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: dHdz
  !internal
  integer :: j
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0.
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
!*  Compute redshift as a function of comoving distance
!*
!*  Args:
!*    :rz[zn] (double) : Comoving distance [Mpc]
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*
!*  Returns:
!*    :z[zn] (double)  : Redshift
!*
  implicit none
  !f2py intent(in) zn, H0, Om, Ov, w0, wa, rz
  !f2py intent(out) z
  !f2py depend(zn) rz, z
  !I/O
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: rz
  double precision, intent(out), dimension(0:zn-1) :: z
  !internal
  integer :: j
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0.
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
!*  Compute comoving distance as a function of z
!*
!*  Args:
!*    :z[zn] (double)  : Redshift
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*
!*  Returns:
!*    :rz[zn] (double) : Comoving distance [Mpc]
!*
  implicit none
  !f2py intent(in) zn, H0, Om, Ov, w0, wa, z
  !f2py intent(out) rz
  !f2py depend(zn) z, rz
  !I/O
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: rz
  !internal
  integer :: j
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0.
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
!*  Compute luminosity distance as a function of z
!*
!*  Args:
!*    :z[zn] (double)  : Redshift
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*
!*  Returns:
!*    :DLz[zn] (double) : Luminosity distance [Mpc]
!*
  implicit none
  !f2py intent(in) zn, H0, Om, Ov, w0, wa, z
  !f2py intent(out) DLz
  !f2py depend(zn) z, DLz
  !I/O
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: DLz
  !internal
  integer :: j
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0.
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
!*  Compute analytic linear growth factor D(z) as a function of z
!*
!*  Args:
!*    :z[zn] (double)  : Redshift
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*    :normed (bool)   : If True, D(z=0)=1. Otherwise, the normalization is defined so that D(z)=a in the pure matter universe, Om(a)=1. 
!*
!*  Returns:
!*    :Dz[zn] (double) : Growth factor
!*
  implicit none
  !f2py intent(in) normed, zn, H0, Om, Ov, w0, wa, z
  !f2py intent(out) Dz
  !f2py depend(zn) z, Dz
  !I/O
  logical, intent(in) :: normed
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: Dz
  !internal
  integer :: j
  double precision :: wz
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0., normed = False
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
!*  Compute linear growth rate f(z) = dlnD/dlna as a function of z
!*
!*  Args:
!*    :z[zn] (double)  : Redshift
!*
!*  Args(optional):
!*    :H0 (double)     : The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
!*    :Om (double)     : The current value of Omega_matter, default to 0.3
!*    :Ov (double)     : The current value of Omega_Dark-energy, default to 0.7
!*    :w0, wa (double) : The EoS of Dark Energy, default to w0=-1 and wa=0.
!*
!*  Returns:
!*    :fz[zn] (double) : Growth rate
!*
  implicit none
  !f2py intent(in) zn, H0, Om, Ov, w0, wa, z
  !f2py intent(out) fz
  !f2py depend(zn) z, fz
  integer, intent(in) :: zn
  double precision, intent(in) :: H0, Om, Ov, w0, wa
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: fz
  integer :: j
  type(cosmoparams) :: cp
  !rmargs :: zn
  !opt4py :: H0 = 70., Ov = 0.7, Om = 0.3, w0 = -1., wa = 0.
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
!*  Distribution function of NS-NS merger events per redshift (dN/dz) at z
!*
!*  Args:
!*    :z (double)  : redshift
!*    :Cz (double) : comoving distance
!*    :Hz (double) : expansion rate
!*
!*  Args(optional):
!*    :ntype (str) : type of dotn functional form, i.e, CH06 (default) or none.
!*    :dotn0 (double) : current merger-rate
!*    :Tobs (double)  : total observation time
!*
!*  Returns:
!*    :nz (double) : distribution function at z
  implicit none
  !f2py intent(in) ntype, z, Cz, Hz, dotn0, Tobs
  !f2py intent(out) nz
  !I/O
  character(4), intent(in) :: ntype
  double precision, intent(in) :: z, Cz, Hz
  double precision, intent(in) :: dotn0, Tobs
  double precision, intent(out) :: nz
  !internal
  double precision :: Lz, sz
  !opt4py :: ntype = 'CH06', dotn0 = 1e-6, Tobs = 3.

  Lz = Cz*(1d0+z)
  sz = dotn(z,ntype)
  nz = Tobs*dotn0*4d0*pi*sz*Lz**2/(Hz*(1d0+z)**3)

end subroutine nz_gw

subroutine drate_dz(z,zn,ntype,dndz)
  implicit none
  !f2py intent(in) ntype, zn, z
  !f2py intent(out) dndz
  !f2py depend(zn) z, dndz
  character(4), intent(in) :: ntype
  integer, intent(in) :: zn
  double precision, intent(in), dimension(0:zn-1) :: z
  double precision, intent(out), dimension(0:zn-1) :: dndz
  integer :: j
  !rmargs :: zn
  !opt4py :: ntype = 'CH06'
  !add2py :: zn = len(z)

  do j = 0, zn-1
    dndz(j) = dlndotn(z(j),ntype)
  end do

end subroutine drate_dz



end module cosmofuncs



