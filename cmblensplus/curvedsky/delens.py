from cmblensplus import libcurvedsky
import numpy

def lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside_t=0,gtype='p'):
  """
  Computing lensing B mode as a convolution of wiener-filtered E-mode and lensing potential

  Args:
    :lmax (*int*): Maximum multipole of output lensing B-mode alm
    :elmin (*int*): Minimum multipole of wiener-filtered E-mode alm
    :elmax (*int*): Maximum multipole of wiener-filtered E-mode alm
    :plmin (*int*): Minimum multipole of wiener-filtered lensing potential alm
    :plmax (*int*): Maximum multipole of wiener-filtered lensing potential alm
    :wElm [*l,m*] (*dcmplx*): Wiener-filtered E-mode alm, with bounds (0:elmax,0:elmax)
    :wplm [*l,m*] (*dcmplx*): Wiener-filtered lensing potential (or kappa) alm, with bounds (0:plmax,0:plmax)

  Args(optional):
    :nside_t (*int*): Nside for the convolution calculation
    :gtype (*str*): Type of input wplm ('p'=phi or 'k'=kappa), default to 'p' (phi).

  Returns:
    :lBlm [*l,m*] (*dcmplx*): Lensing B-mode alm, with bounds (0:lmax,0:lmax)

  Usage:
    :lBlm = curvedsky.delens.lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside_t,gtype):
  """
  return libcurvedsky.delens.lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside_t,gtype)

def shiftvec(nside,lmax,plm,nremap):
  """
  Return the anti deflection vector, beta, at the Healpix pixel for the delensing where 

    beta(n) + alphaiw(n+beta(n)) = 0

  and alphaw is the filtered lensing deflection vector (see arXiv:1701.01712).

  Args:
    :nside (*int*): Nside of output shift vector
    :lmax (*int*): Maximum multipole of the input plm
    :plm [*l,m*] (*dcmplx*): Wiener-filtered lensing potential alm, with bounds (0:lmax,0:lmax)

  Args(optional):
    :nremap (*int*): Number of iteration for computing the shift vector

  Returns:
    :beta [*pix,2*] (*double*): 2D shift vector, with bounds (0:npix-1,1:2)

  Usage:
    :beta = curvedsky.delens.shiftvec(nside,lmax,plm,nremap):
  """
  npix = 12*nside**2
  return libcurvedsky.delens.shiftvec(npix,lmax,plm,nremap)

def phi2grad(nside,lmax,plm):
  """
  Return the deflection vector, grad, at the Healpix pixel

  Args:
    :nside (*int*): Nside of output deflection vector
    :lmax (*int*): Maximum multipole of the input plm/clm
    :plm [*l,m*] (*dcmplx*): Lensing potential alm, with bounds (0:lmax,0:lmax)

  Returns:
    :grad [*pix,2*] (*double*): 2D deflection vector, with bounds (0:npix-1,1:2)

  Usage:
    :grad = curvedsky.delens.phi2grad(nside,lmax,plm):
  """
  npix = 12*nside**2
  return libcurvedsky.delens.phi2grad(npix,lmax,plm)

def remap_tp(nside,lmax,beta,alm_in):
  """
  Remapping CMB temperaure and polarization with a given shift vector based on a simple implementation of LensPix
  This function returs X(n+beta) where n is the original direction and beta is the shift vector
  The output is given by alms where alm[0,l,m] is temperature, alm[1,l,m] is E mode, and alm[2,l,m] is B mode.

  Args:
    :nside (*int*): Nside of input shift vector
    :lmax (*int*): Maximum multipole of the input plm
    :beta [*pix,2*] (*double*): 2D shift vector, with bounds (0:npix-1,1:2)
    :alm_in [*TEB,l,m*] (*dcmplx*): Input T/E/B alms to be remapped, with bounds (1:3,0:lmax,0:lmax).

  Returns:
    :alm_re [*TEB,l,m*] (*dcmplx*): Remapped T/E/B alms, with bounds (1:3,0:lmax,0:lmax).

  Usage:
    :alm_re = curvedsky.delens.remap_tp(nside,lmax,beta,alm_in):
  """
  npix = 12*nside**2
  return libcurvedsky.delens.remap_tp(npix,lmax,beta,alm_in)

def remap_tp_map(nside,lmax,beta,alm_in):
  """
  Remapping CMB temperaure and polarization with a given shift vector based on a simple implementation of LensPix
  This function returs X(n+beta) where n is the original direction and beta is the shift vector
  The output is given by alms where alm[0,l,m] is temperature, alm[1,l,m] is E mode, and alm[2,l,m] is B mode.

  Args:
    :nside (*int*): Nside of input shift vector
    :lmax (*int*): Maximum multipole of the input plm
    :beta [*pix,2*] (*double*): 2D shift vector, with bounds (0:npix-1,1:2)
    :alm_in [*TEB,l,m*] (*dcmplx*): Input T/E/B alms to be remapped, with bounds (1:3,0:lmax,0:lmax).

  Returns:
    :map_re [*pix,TQU*] (*double*): Remapped T/Q/U alms, with bounds (0:lmax,0:lmax,1:3).

  Usage:
    :map_re = curvedsky.delens.remap_tp_map(nside,lmax,beta,alm_in):
  """
  npix = 12*nside**2
  return libcurvedsky.delens.remap_tp_map(npix,lmax,beta,alm_in)

