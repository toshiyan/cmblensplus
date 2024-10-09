import libbasic
import numpy

def alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,gln= 100,lxcut= 0,gle= 1e-14):
  """
  Compute flat-sky quadratic estimator normalization
  CAUTION: This code interpolates the input Cl at the non-integer multipole by simply Cl(int(ell)) which leads to a small discrepancy in the normalization computed from the FFT-based method (which uses linear interpolation) and from this code. It is desireble to use the FFT-based normalization if you want to normalize the simulation results.
 
  Args:
    :qest (*str*): estimator combination (TT, TE, TB, EE, EB, or BB)
    :qtype (*str*): estimator type (lensing, patchytau)
    :lmax (*double*): output maximum multipole
    :rlmax/rlmin (*double*): input CMB multipole range for reconstruction
    :fC[*rlmax*] (*double*): power spectrum in the numerator
    :W1/W2[*rlmax*] : inverse of the observed power spectrum

  Args(optional):
    :gln (*int*): number of the GL integration points
    :lxcut (*int*): multipole cut in x-direction, |l_x| < lx
    :gle (*double*): convergence parameter for the GL integration

  Returns:
    :Ag/Ac[*l*] (*double*): normalization for even and odd estimator pairs

  Usage:
    :Ag,Ac = basic.flat.alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,gln,gle,lxcut):
  """
  return libbasic.flat.alxy(qest,qtype,lmax,rlmin,rlmax,fC,W1,W2,gln,gle,lxcut)

def alxy_asym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,gln= 100,lxcut= 0,gle= 1e-14):
  """
  Usage:
    :Ag,Ac = basic.flat.alxy_asym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,gln,gle,lxcut):
  """
  return libbasic.flat.alxy_asym(qest,qtype,lmax,rlmin,rlmax,fC,AA,BB,AB,gln,gle,lxcut)

def bbxy(lmax,rlmin,rlmax,XX,YY,weight= 'lensing',gln= 100,gle= 1e-14):
  """
  Usage:
    :BB = basic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle):
  """
  return libbasic.flat.bbxy(lmax,rlmin,rlmax,XX,YY,weight,gln,gle)

