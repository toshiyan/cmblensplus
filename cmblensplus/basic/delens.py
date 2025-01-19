from cmblensplus import libbasic
import numpy

def lintemplate(lmax,elmin,elmax,klmin,klmax,CE,Cm,WE,Wm,gtype='p'):
  """
  Estimate of lensing template B-mode power spectrum (Wiener filters as inputs)

  Args:
    :lmax (*int*): Maximum multipole of output spectrum
    :elmin (*int*): Minimum multipole of E
    :elmax (*int*): Maximum multipole of E
    :klmin (*int*): Minimum multipole of lensing mass
    :klmax (*int*): Maximum multipole of lensing mass
    :CE[*l*] (*double*): Power spectrum of E-mode, with bounds (0:dlmax)
    :Cp[*l*] (*double*): Power spectrum of lensing pontential, with bounds (0:dlmax)
    :WE[*l*] (*double*): Wiener filter of E-mode, with bounds (0:dlmax)
    :Wp[*l*] (*double*): Wiener filter of lensing potential, with bountd (0:dlmax)

  Args(optional):
    :gtype (*str*): specify type of the input Cp, p (default) or k.

  Returns:
    :CB[*l*] (*double*): Lensing B-mode power spectrum, with bounds (0:lmax)

  Usage:
    :CB = basic.delens.lintemplate(lmax,elmin,elmax,klmin,klmax,CE,Cm,WE,Wm,gtype):
  """
  return libbasic.delens.lintemplate(lmax,elmin,elmax,klmin,klmax,CE,Cm,WE,Wm,gtype)

def lensingbb(lmax,dlmin,dlmax,CE,Cp):
  """
 Lensing B-mode power spectrum as a convolution of ClEE and Clpp

  Args:
    :lmax (*int*): Maximum multipole of output spectrum
    :dlmin (*int*): Minimum multipole of E and lensing for delensing
    :dlmax (*int*): Maximum multipole of E and lensing for delensing
    :CE[*l*] (*double*): Power spectrum of E-mode, with bounds (0:dlmax)
    :Cp[*l*] (*double*): Power spectrum of lensing pontential, with bounds (0:dlmax)

  Returns:
    :CB[*l*] (*double*): Lensing B-mode power spectrum, with bounds (0:lmax)

  Usage:
    :CB = basic.delens.lensingbb(lmax,dlmin,dlmax,CE,Cp):
  """
  return libbasic.delens.lensingbb(lmax,dlmin,dlmax,CE,Cp)

def delensbias_dom(lmax,dlmin,dlmax,CE,CB,Cp,NP1,NP2,Ag):
  """
  Dominant term of the delensing bias in the B-mode internal delensing

  Args:
    :lmax (*int*): Maximum multipole of output spectrum
    :dlmin (*int*): Minimum multipole of E and lensing for delensing
    :dlmax (*int*): Maximum multipole of E and lensing for delensing
    :CE[*l*] (*double*): Power spectrum of E-mode, with bounds (0:dlmax)
    :CB[*l*] (*double*): Power spectrum of B-mode, with bounds (0:dlmax)
    :Cp[*l*] (*double*): Power spectrum of lensing pontential, with bounds (0:dlmax)
    :NP1[*l*] (*double*): Pol. noise spectrum for lensing reconstruction, with bounds (0:dlmax)
    :NP2[*l*] (*double*): Pol. noise spectrum for B-mode to be delensed, with bounds (0:dlmax)
    :Ag[*l*] (*double*): Lensing reconstruction noise, with bounds (0:dlmax)

  Returns:
    :DB[*l*] (*double*): Lensing B-mode power spectrum, with bounds (0:lmax)

  Usage:
    :DB = basic.delens.delensbias_dom(lmax,dlmin,dlmax,CE,CB,Cp,NP1,NP2,Ag):
  """
  return libbasic.delens.delensbias_dom(lmax,dlmin,dlmax,CE,CB,Cp,NP1,NP2,Ag)

