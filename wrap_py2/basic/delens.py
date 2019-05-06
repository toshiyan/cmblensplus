import libbasic

def resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp):
  """
  Residual B-mode spectrum; ClBB = ClBB^lin - ClBB^est

  Args:
    lmax (int)     : maximum multipole of residual ClBB
    dlmin (int)    : minimum multipole of E and lensing for delensing
    dlmax (int)    : maximum multipole of E and lensing for delensing
    CE[l] (double) : power spectrum of E-mode, with bounds (0:dlmax)
    Cp[l] (double) : power spectrum of lensing pontential, with bounds (0:dlmax)
    WE[l] (double) : Wiener filter of E-mode, with bounds (0:dlmax)
    Wp[l] (double) : Wiener filter of lensing potential, with bountd (0:dlmax)

  Returns:
    CB[l] (double) : residual B-mode spectrum, with bounds (0:lmax)

  Usage:
    - e.g., lmax,CB = basic.delens.resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp)
  """
  return libbasic.delens.resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp)

def lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp):
  """
  Usage:
    - e.g., lmax,Cl = basic.delens.lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp)
  """
  return libbasic.delens.lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp)

def lensingbb(lmax,dlmin,dlmax,CE,Cp):
  """
  Usage:
    - e.g., lmax,CB = basic.delens.lensingbb(lmax,dlmin,dlmax,CE,Cp)
  """
  return libbasic.delens.lensingbb(lmax,dlmin,dlmax,CE,Cp)

def delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl):
  """
  Usage:
    - e.g., lmax = basic.delens.delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)
  """
  return libbasic.delens.delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)

