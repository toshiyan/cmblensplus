import libanalytic

def qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter=0,conv=0):
  """
     lensing reconstruction with EB
     convergence check
     delensing with EB-estimator
  Usage:
    - e.g., Alg,Alc = analytic.delens.qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter,conv)
  """
  if iter==0: iter= 1
  if conv==0: conv= 0.001
  return analytic.delens.qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,iter,conv)

def resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp):
  """
 residual ClBB = ClBB^lin - ClBB^est
  Usage:
    - e.g., lmax,CB = analytic.delens.resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp)
  """
  return analytic.delens.resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp)

def lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp):
  """
  Usage:
    - e.g., lmax,Cl = analytic.delens.lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp)
  """
  return analytic.delens.lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp)

def lensingbb(lmax,dlmin,dlmax,CE,Cp):
  """
  Usage:
    - e.g., lmax,CB = analytic.delens.lensingbb(lmax,dlmin,dlmax,CE,Cp)
  """
  return analytic.delens.lensingbb(lmax,dlmin,dlmax,CE,Cp)

def delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl):
  """
  Usage:
    - e.g., lmax = analytic.delens.delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)
  """
  return analytic.delens.delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)

