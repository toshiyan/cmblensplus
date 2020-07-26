import libcurvedsky

def qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t=0,verbose=False):
  """
  Reconstructing inhomogeneous tau from the temperature quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output tau alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fC [*l*] (*double*): TT spectrum, with bounds (1:rlmax)
    :Tlm1 [*l,m*] (*dcmplx*): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
    :Tlm2 [*l,m*] (*dcmplx*): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :verbose (*bool*): Output messages, default to False

  Returns:
    :alm [*l,m*] (*dcmplx*): Amplitude modulation alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.rec_tau.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t,verbose):
  """
  return libcurvedsky.rec_tau.qtt(lmax,rlmin,rlmax,fC,Tlm1,Tlm2,nside_t,verbose)

def qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,nside_t=0,verbose=False):
  """
  Reconstructing amplitude modulation from the EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fCE [*l*] (*double*): EE spectrum, with bounds (0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :verbose (*bool*): Output messages, default to False

  Returns:
    :alm [*l,m*] (*dcmplx*): Amplitude modulation alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.rec_tau.qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,nside_t,verbose):
  """
  return libcurvedsky.rec_tau.qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,nside_t,verbose)

def oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,nside_t=0,verbose=False):
  """
  Reconstructing amplitude modulation by the odd EB quadratic estimator

  Args:
    :lmax (*int*): Maximum multipole of output lensing potential alms
    :rlmin/rlmax (*int*): Minimum/Maximum multipole of CMB for reconstruction
    :fEB [*l*] (*double*): EB spectrum, with bounds (0:rlmax)
    :Elm [*l,m*] (*dcmplx*): Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
    :Blm [*l,m*] (*dcmplx*): Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :verbose (*bool*): Output messages, default to False

  Returns:
    :alm [*l,m*] (*dcmplx*): Reconstructed alm, with bounds (0:lmax,0:lmax)

  Usage:
    :alm = curvedsky.rec_tau.oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,nside_t,verbose):
  """
  return libcurvedsky.rec_tau.oeb(lmax,rlmin,rlmax,fEB,Elm,Blm,nside_t,verbose)

