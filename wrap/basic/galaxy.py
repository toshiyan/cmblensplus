import libbasic

def zbin(zn,a,b,zm,verbose=False):
  """
 Computing z-interval of z-bin so that number of galaxies at each z-bin is equal

  Args:
    :zn (*int*): number of z-bins
    :a, b (*double*): shape parameters of Schechter-like galaxy distribution
    :zm (*double*): mean redshift

  Args(optional):
    :verbose (*bool*): output messages (default to True)

  Returns:
    :zb[*zn+1*] (*double*): z-intervals

  Usage:
    :zb = basic.galaxy.zbin(zn,a,b,zm,verbose):
  """
  return libbasic.galaxy.zbin(zn,a,b,zm,verbose)

def frac(zn,zb,a,b,zm,zbias=0.0,sigma=0.0,verbose=False):
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

