import libbasic

def dndz_sf(z,a,b,zm=0,z0=0,zn=None):
  """
 A model of galaxy z distribution

  Args:
    :z[*zn*] (*double*): redshifts at which dNdz is returned
    :a, b (*double*): shape parameters of Schechter-like galaxy distribution

  Args(optional):
    :zm (*double*): mean redshift, default to 0
    :z0 (*double*): a parameter relevant to zm, default to 0. Either zm or z0 has to be specified.

  Returns:
    :dndz[*zn*] (*double*): galaxy z distribution

  Usage:
    :dndz = basic.galaxy.dndz_sf(zn,z,a,b,zm,z0):
  """
  if zn is None: zn=len(z)
  return libbasic.galaxy.dndz_sf(zn,z,a,b,zm,z0)

def photoz_error(z,zi,zn=None,sigma=0.03,zbias=0.):
  """
 Photo-z error on z distribution which is multiplied to original galaxy z distribution. 
 See Eq.(13) of arXiv:1103.1118 for details.

  Args:
    :z[*zn*] (*double*): redshifts at which photoz error function is returned
    :zi[*2*] (*double*): z-bin edges
    :sigma (*double*): a parameter of photo-z error which is given by, sigma x (1+z)
    :zbias (*double*): photo-z mean bias

  Returns:
    :pz[*zn*] (*double*): photoz error function

  Usage:
    :pz = basic.galaxy.photoz_error(zn,z,zi,sigma,zbias):
  """
  if zn is None: zn=len(z)
  return libbasic.galaxy.photoz_error(zn,z,zi,sigma,zbias)

def zbin(zn,a,b,zm=0,z0=0,verbose=False):
  """
 Computing z-interval of z-bin so that number of galaxies at each z-bin is equal

  Args:
    :zn (*int*): number of z-bins
    :a, b (*double*): shape parameters of Schechter-like galaxy distribution

  Args(optional):
    :zm (*double*): mean redshift, default to 0
    :z0 (*double*): a parameter relevant to zm, default to 0. Either zm or z0 has to be specified.
    :verbose (*bool*): output messages (default to True)

  Returns:
    :zb[*zn+1*] (*double*): z-intervals

  Usage:
    :zb = basic.galaxy.zbin(zn,a,b,zm,z0,verbose):
  """
  return libbasic.galaxy.zbin(zn,a,b,zm,z0,verbose)

def frac(zn,zb,a,b,zm,verbose=False,zbias=0.0,sigma=0.0):
  """
 Computing z-interval of z-bin so that number of galaxies at each z-bin is equal

  Args:
    :zn (*int*): number of z-bins
    :a, b (*double*): shape parameters of Schechter-like galaxy distribution
    :zm (*double*): mean redshift
    :zb[*zn+1*] (*double*): z-intervals

  Args(optional):
    :verbose (*bool*): output messages (default to True)
    :zbias (*double*): constant bias to true photo-z
    :sigma (*double*): uncertaines of photo-z

  Returns:
    :nfrac[*zn*] (*double*): fraction of galaxy number at each bin, defined by int_zi^zi+1 dz N(z)/int dz N(z)

  Usage:
    :nfrac = basic.galaxy.frac(zn,zb,a,b,zm,zbias,sigma,verbose):
  """
  return libbasic.galaxy.frac(zn,zb,a,b,zm,zbias,sigma,verbose)

