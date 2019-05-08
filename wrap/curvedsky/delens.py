import libcurvedsky

def lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside=None):
  """
  Computing lensing B mode as a convolution of wiener-filtered E-mode and lensing potential

  Args:
    :lmax (*int*): Maximum multipole of output lensing B-mode alm
    :elmin (*int*): Minimum multipole of wiener-filtered E-mode alm
    :elmax (*int*): Maximum multipole of wiener-filtered E-mode alm
    :plmin (*int*): Minimum multipole of wiener-filtered lensing potential alm
    :plmax (*int*): Maximum multipole of wiener-filtered lensing potential alm
    :wElm [*l,m*] (*dcmplx*): Wiener-filtered E-mode alm, with bounds (0:elmax,0:elmax)
    :wplm [*l,m*] (*dcmplx*): Wiener-filtered lensing potential alm, with bounds (0:plmax,0:plmax)

  Args(optional):
    :nside (*int*): Nside for the convolution calculation, default to lmax

  Returns:
    :lBlm [*l,m*] (*dcmplx*): Lensing B-mode alm, with bounds (0:lmax,0:lmax)

  Usage:
    :lBlm = curvedsky.delens.lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside):
  """
  if nside is None: nside= lmax
  return libcurvedsky.delens.lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside)

def shiftvec(npix,lmax,plm,nremap= 3):
  """
  Return the anti deflection vector, beta, at the Healpix pixel for the delensing where 

    beta(n) + alphaiw(n+beta(n)) = 0

  and alphaw is the filtered lensing deflection vector (see arXiv:1701.01712).

  Args:
    :npix (*int*): Pixel number of output shift vector
    :lmax (*int*): Maximum multipole of the input plm
    :plm [*l,m*] (*dcmplx*): Wiener-filtered lensing potential alm, with bounds (0:lmax,0:lmax)

  Args(optional):
    :nremap (*int*): Number of iteration for computing the shift vector

  Returns:
    :beta [*pix,2*] (*double*): 2D shift vector, with bounds (0:npix-1,1:2)

  Usage:
    :beta = curvedsky.delens.shiftvec(npix,lmax,plm,nremap):
  """
  return libcurvedsky.delens.shiftvec(npix,lmax,plm,nremap)

def phi2grad(npix,lmax,plm):
  """
  Return the deflection vector, grad, at the Healpix pixel

  Args:
    :npix (*int*): Pixel number of output deflection vector
    :lmax (*int*): Maximum multipole of the input plm/clm
    :plm [*l,m*] (*dcmplx*): Lensing potential alm, with bounds (0:lmax,0:lmax)

  Returns:
    :grad [*pix,2*] (*double*): 2D deflection vector, with bounds (0:npix-1,1:2)

  Usage:
    :grad = curvedsky.delens.phi2grad(npix,lmax,plm):
  """
  return libcurvedsky.delens.phi2grad(npix,lmax,plm)

def remap_tp(npix,lmax,beta,alm_in):
  """
  Remapping CMB temperaure and polarization with a given shift vector, grad, based on a simple implementation of LensPix. 
  alm[*0,:,:*] is temperature, alm[*1,:,:*] is E mode, and alm[*2,:,:*] is B mode.

  Args:
    :npix (*int*): Pixel number of output shift vector
    :lmax (*int*): Maximum multipole of the input plm
    :beta [*pix,2*] (*double*): 2D shift vector, with bounds (0:npix-1,1:2)
    :alm_in [*TEB,l,m*] (*dcmplx*): Input T/E/B alms to be remapped, with bounds (1:3,0:lmax,0:lmax).

  Returns:
    :alm_re [*TEB,l,m*] (*dcmplx*): Remapped T/E/B alms, with bounds (1:3,0:lmax,0:lmax).

  Usage:
    :alm_re = curvedsky.delens.remap_tp(npix,lmax,beta,alm_in):
  """
  return libcurvedsky.delens.remap_tp(npix,lmax,beta,alm_in)

