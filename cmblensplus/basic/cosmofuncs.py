from cmblensplus import libbasic
import numpy

def hubble(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,divc=False):
  """
  Compute the expansion rate in unit of 1/Mpc, H/c, or in unit of km/s/Mpc, H.  

  Args:
    :z[*zn*] (*double*): Redshifts at which H is computed

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.
    :divc (*bool*): Divide H by c or not, default to False.

  Returns:
    :Hz[*zn*] (*double*): The expansion rate, H(z), divided by c or not.

  Usage:
    :Hz = basic.cosmofuncs.hubble(z,H0,Om,Ov,w0,wa,zn,divc):
  """
  zn = len(z)
  return libbasic.cosmofuncs.hubble(z,H0,Om,Ov,w0,wa,zn,divc)

def dhubble_dz(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.):
  """
  Compute dH(z)/dz.  

  Args:
    :z[*zn*] (*double*): Redshifts at which dH/dz is computed

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.

  Returns:
    :dHdz[*zn*] (*double*): The derivative of the expansion rate, dH(z)/dz.

  Usage:
    :dHdz = basic.cosmofuncs.dhubble_dz(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.dhubble_dz(z,H0,Om,Ov,w0,wa,zn)

def dist2z(rz,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.):
  """
  Compute redshift as a function of comoving distance

  Args:
    :rz[*zn*] (*double*): Comoving distance [*Mpc*]

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.

  Returns:
    :z[*zn*] (*double*): Redshift

  Usage:
    :z = basic.cosmofuncs.dist2z(rz,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(rz)
  return libbasic.cosmofuncs.dist2z(rz,H0,Om,Ov,w0,wa,zn)

def dist_comoving(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.):
  """
  Compute comoving distance as a function of z

  Args:
    :z[*zn*] (*double*): Redshift

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.

  Returns:
    :rz[*zn*] (*double*): Comoving distance [*Mpc*]

  Usage:
    :rz = basic.cosmofuncs.dist_comoving(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.dist_comoving(z,H0,Om,Ov,w0,wa,zn)

def dist_luminosity(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.):
  """
  Compute luminosity distance as a function of z

  Args:
    :z[*zn*] (*double*): Redshift

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.

  Returns:
    :DLz[*zn*] (*double*): Luminosity distance [*Mpc*]

  Usage:
    :DLz = basic.cosmofuncs.dist_luminosity(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.dist_luminosity(z,H0,Om,Ov,w0,wa,zn)

def growth_factor(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,normed=False):
  """
  Compute analytic linear growth factor D(z) as a function of z

  Args:
    :z[*zn*] (*double*): Redshift

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.
    :normed (*bool*): If True, D(z=0)=1. Otherwise, the normalization is defined so that D(z)=a in the pure matter universe, Om(a)=1.

  Returns:
    :Dz[*zn*] (*double*): Growth factor

  Usage:
    :Dz = basic.cosmofuncs.growth_factor(z,H0,Om,Ov,w0,wa,zn,normed):
  """
  zn = len(z)
  return libbasic.cosmofuncs.growth_factor(z,H0,Om,Ov,w0,wa,zn,normed)

def growth_rate(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.):
  """
  Compute linear growth rate f(z) = dlnD/dlna as a function of z

  Args:
    :z[*zn*] (*double*): Redshift

  Args(optional):
    :H0 (*double*): The current value of hubble parameter in km/s/Mpc, default to 70 km/s/Mpc
    :Om (*double*): The current value of Omega_matter, default to 0.3
    :Ov (*double*): The current value of Omega_Dark-energy, default to 0.7
    :w0, wa (*double*): The EoS of Dark Energy, default to w0=-1 and wa=0.

  Returns:
    :fz[*zn*] (*double*): Growth rate

  Usage:
    :fz = basic.cosmofuncs.growth_rate(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.growth_rate(z,H0,Om,Ov,w0,wa,zn)

def nz_gw(z,Cz,Hz,ntype='CH06',dotn0=1e-6,Tobs=3.):
  """
  Distribution function of NS-NS merger events per redshift (dN/dz) at z

  Args:
    :z (*double*): redshift
    :Cz (*double*): comoving distance
    :Hz (*double*): expansion rate

  Args(optional):
    :ntype (*str*): type of dotn functional form, i.e, CH06 (default) or none.
    :dotn0 (*double*): current merger-rate
    :Tobs (*double*): total observation time

  Returns:
    :nz (*double*): distribution function at z
  Usage:
    :nz = basic.cosmofuncs.nz_gw(z,Cz,Hz,ntype,dotn0,Tobs):
  """
  return libbasic.cosmofuncs.nz_gw(z,Cz,Hz,ntype,dotn0,Tobs)

def drate_dz(z,ntype='CH06'):
  """
  Usage:
    :dndz = basic.cosmofuncs.drate_dz(z,zn,ntype):
  """
  zn = len(z)
  return libbasic.cosmofuncs.drate_dz(z,zn,ntype)

