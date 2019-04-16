import curvedsky

def lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside=0):
  """
  Computing lensing B mode as a convolution of wiener-filtered E-mode and lensing potential

  Args:
    - lmax (int)        : maximum multipole of output lensing B-mode alm
    - elmin (int)       : minimum multipole of wiener-filtered E-mode alm
    - elmax (int)       : maximum multipole of wiener-filtered E-mode alm
    - plmin (int)       : minimum multipole of wiener-filtered lensing potential alm
    - plmax (int)       : maximum multipole of wiener-filtered lensing potential alm
    - wElm[l,m] (dcmplx): wiener-filtered E-mode alm, with bounds (0:elmax,0:elmax)
    - wplm[l,m] (dcmplx): wiener-filtered lensing potential alm, with bounds (0:plmax,0:plmax)

  Args(optional):
    - nside (int)       : Nside for the convolution calculation, default to lmax

  Returns:
    - lBlm[l,m] (dcmplx): Lensing B-mode alm, with bounds (0:lmax,0:lmax)

  Usage:
    - e.g., lBlm = curvedsky.delens.lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside)
  """
  if nside==0: nside= lmax
  return curvedsky.delens.lensingb(lmax,elmin,elmax,plmin,plmax,wElm,wplm,nside)

def shiftvec(npix,lmax,plm,nremap=0):
  """
  Return the anti deflection vector, beta, at the Healpix pixel for the delensing where 

    beta(n) + alphaiw(n+beta(n)) = 0

  and alphaw is the filtered lensing deflection vector (see arXiv:1701.01712).

  Args:
    - npix (int)          : pixel number of output shift vector
    - lmax (int)          : maximum multipole of the input plm
    - plm[l,m] (dcmplx)   : wiener-filtered lensing potential alm, with bounds (0:lmax,0:lmax)

  Args(optional):
    - nremap (int)        : number of iteration for computing the shift vector

  Returns:
    - beta[pix,2] (double): 2D shift vector, with bounds (0:npix-1,1:2)

  Usage:
    - e.g., beta = curvedsky.delens.shiftvec(npix,lmax,plm,nremap)
  """
  if nremap==0: nremap= 3
  return curvedsky.delens.shiftvec(npix,lmax,plm,nremap)

def remap_tp(npix,lmax,beta,alm_in):
  """
  Remapping CMB temperaure and polarization with a given shift vector, grad, based on a simple implementation of LensPix. 
  alm[0,:,:] is temperature, alm[1,:,:] is E mode, and alm[2,:,:] is B mode.

  Args:
    - npix (int)               : pixel number of output shift vector
    - lmax (int)               : maximum multipole of the input plm
    - beta[pix,2] (double)     : 2D shift vector, with bounds (0:npix-1,1:2)
    - alm_in[TEB,l,m] (dcmplx) : input T/E/B alms to be remapped, with bounds (1:3,0:lmax,0:lmax).

  Returns:
    - alm_re[TEB,l,m] (dcmplx) : remapped T/E/B alms, with bounds (1:3,0:lmax,0:lmax). 

  Usage:
    - e.g., alm_re = curvedsky.delens.remap_tp(npix,lmax,beta,alm_in)
  """
  return curvedsky.delens.remap_tp(npix,lmax,beta,alm_in)

