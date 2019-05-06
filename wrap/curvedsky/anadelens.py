import libcurvedsky

def qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter=None,conv=None):
  """
     lensing reconstruction with EB
     convergence check
     delensing with EB-estimator
  Usage:
    - e.g., Alg,Alc = curvedsky.anadelens.qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter,conv)
  """
  if iter is None: iter= 1
  if conv is None: conv= 0.001
  return libcurvedsky.anadelens.qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter,conv)

def resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp):
  """
 residual ClBB = ClBB^lin - ClBB^est
  Usage:
    - e.g., lmax,CB = curvedsky.anadelens.resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp)
  """
  return libcurvedsky.anadelens.resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp)

def lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp):
  """
  Usage:
    - e.g., lmax,Cl = curvedsky.anadelens.lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp)
  """
  return libcurvedsky.anadelens.lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp)

def lensingbb(lmax,dlmin,dlmax,CE,Cp):
  """
  Usage:
    - e.g., lmax,CB = curvedsky.anadelens.lensingbb(lmax,dlmin,dlmax,CE,Cp)
  """
  return libcurvedsky.anadelens.lensingbb(lmax,dlmin,dlmax,CE,Cp)

def delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl):
  """
  Usage:
    - e.g., lmax = curvedsky.anadelens.delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)
  """
  return libcurvedsky.anadelens.delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)

