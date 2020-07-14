import libbasic

def hubble(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0):
  """
  Usage:
    :Hz = basic.cosmofuncs.hubble(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.hubble(z,H0,Om,Ov,w0,wa,zn)

def dhubble_dz(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0):
  """
  Usage:
    :dHdz = basic.cosmofuncs.dhubble_dz(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.dhubble_dz(z,H0,Om,Ov,w0,wa,zn)

def dist2z(rz,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0):
  """
  Usage:
    :z = basic.cosmofuncs.dist2z(rz,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(rz)
  return libbasic.cosmofuncs.dist2z(rz,H0,Om,Ov,w0,wa,zn)

def dist_comoving(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0):
  """
  Usage:
    :rz = basic.cosmofuncs.dist_comoving(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.dist_comoving(z,H0,Om,Ov,w0,wa,zn)

def dist_luminosity(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0):
  """
  Usage:
    :DLz = basic.cosmofuncs.dist_luminosity(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.dist_luminosity(z,H0,Om,Ov,w0,wa,zn)

def growth_factor(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0,normed=False):
  """
  Usage:
    :Dz = basic.cosmofuncs.growth_factor(z,H0,Om,Ov,w0,wa,zn,normed):
  """
  zn = len(z)
  return libbasic.cosmofuncs.growth_factor(z,H0,Om,Ov,w0,wa,zn,normed)

def growth_rate(z,H0=70.,Om=0.3,Ov=0.7,w0=-1.,wa=0.,zn=0):
  """
  Usage:
    :fz = basic.cosmofuncs.growth_rate(z,H0,Om,Ov,w0,wa,zn):
  """
  zn = len(z)
  return libbasic.cosmofuncs.growth_rate(z,H0,Om,Ov,w0,wa,zn)

def nz_gw(z,Cz,Hz,ntype,dotn0=1e-6,Tobs=3.):
  """
  Distribution function of NS binary sources per redshift (dN/dz)
  Usage:
    :nz = basic.cosmofuncs.nz_gw(z,Cz,Hz,ntype,dotn0,Tobs):
  """
  return libbasic.cosmofuncs.nz_gw(z,Cz,Hz,ntype,dotn0,Tobs)

def drate_dz(z,zn=0,ntype='CH06'):
  """
  Usage:
    :dndz = basic.cosmofuncs.drate_dz(z,zn,ntype):
  """
  zn = len(z)
  return libbasic.cosmofuncs.drate_dz(z,zn,ntype)

