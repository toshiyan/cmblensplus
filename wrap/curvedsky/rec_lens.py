import libcurvedsky
import numpy

def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): TT spectrum, with bounds (0:rlmax)
    :Tlm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Tlm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (*bool*): Output messages, default to False

  Returns:
    :glm [*l,m*] (*dcmplx*): CMB lensing potential alm, with bounds (0:lmax,0:lmax)
    :clm [*l,m*] (*dcmplx*): Curl mode (pseudo lensing potential) alm, with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t,gtype,verbose):
  """
  return libcurvedsky.rec_lens.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t,gtype,verbose)

def qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the TE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): TE spectrum, with bounds (0:rlmax)
    :Tlm [*l,m*] (*dcmplx*): Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (*bool*): Output messages, default to False

  Returns:
    :glm [*l,m*] (*dcmplx*): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [*l,m*] (*dcmplx*): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.rec_lens.qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside_t,gtype,verbose):
  """
  return libcurvedsky.rec_lens.qte(lmax,rlmin,rlmax,fC,Tlm,Elm,nside_t,gtype,verbose)

def qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the TB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): TE spectrum, with bounds (0:rlmax)
    :Tlm [*l,m*] (*dcmplx*): Inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (*bool*): Output messages, default to False

  Returns:
    :glm [*l,m*] (*dcmplx*): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [*l,m*] (*dcmplx*): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.rec_lens.qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside_t,gtype,verbose):
  """
  return libcurvedsky.rec_lens.qtb(lmax,rlmin,rlmax,fC,Tlm,Blm,nside_t,gtype,verbose)

def qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the EE quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): EE spectrum, with bounds (0:rlmax)
    :Elm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Elm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (*bool*): Output messages, default to False

  Returns:
    :glm [*l,m*] (*dcmplx*): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [*l,m*] (*dcmplx*): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.rec_lens.qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside_t,gtype,verbose):
  """
  return libcurvedsky.rec_lens.qee(lmax,rlmin,rlmax,fC,Elm1,Elm2,nside_t,gtype,verbose)

def qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): EE spectrum, with bounds (0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (*bool*): Output messages, default to False

  Returns:
    :glm [*l,m*] (*dcmplx*): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [*l,m*] (*dcmplx*): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.rec_lens.qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside_t,gtype,verbose):
  """
  return libcurvedsky.rec_lens.qeb(lmax,rlmin,rlmax,fC,Elm,Blm,nside_t,gtype,verbose)

def qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside_t=0,gtype='',verbose=False):
  """
  Reconstructing CMB lensing potential and its curl mode from the BB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipoles of CMB for reconstruction
    :fC [*l*] (*double*): BB spectrum, with bounds (0:rlmax)
    :Blm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of output, i.e., convergence (gtype='k') or lensing potential (gtype='', default)
    :verbose (*bool*): Output messages, default to False

  Returns:
    :glm [*l,m*] (*dcmplx*): CMB lensing potential, with bounds (0:lmax,0:lmax)
    :clm [*l,m*] (*dcmplx*): Curl mode (pseudo lensing potential), with bounds (0:lmax,0:lmax)

  Usage:
    :glm,clm = curvedsky.rec_lens.qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside_t,gtype,verbose):
  """
  return libcurvedsky.rec_lens.qbb(lmax,rlmin,rlmax,fC,Blm1,Blm2,nside_t,gtype,verbose)

