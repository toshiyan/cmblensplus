import libcurvedsky

def cnfilter(n,npix,lmax,cl,nij,alm,itern,eps=1e-6,filter=''):
  """
 Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
 This code assumes
   - The signal power spectrum is isotropic Gaussian. 
   - Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
   - The data model is bxS+N

  Args:
    :n (*int*): Number of maps
    :npix (*int*): Number of pixel
    :lmax (*int*): Maximum multipole of alm
    :cl[*n,l*] (*double*): Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
    :nij[*n,pix*] (*double*): Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
    :alm[*n,l,m*] (*dcmplx*): Input alm, with bouds (0:n-1,0:lmax,0:lmax)
    :itern (*int*): Number of interation

  Args(optional): 
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
    :filter (*str*): C-inverse ('') or Wiener filter (W), default to C-inverse.

  Returns:
    :xlm[*n,l,m*] (*dcmplx*): C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.cninv.cnfilter(n,npix,lmax,cl,nij,alm,itern,eps,filter):
  """
  return libcurvedsky.cninv.cnfilter(n,npix,lmax,cl,nij,alm,itern,eps,filter)

def cg_algorithm(n,npix,lmax,clh,nij,b,itern,eps=1e-6):
  """
  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
  The code assumes 
    - A = [1 + C^1/2 N^-1 C^1/2]
    - C^1/2 is diagonal
    - N is diagonal in pixel space (statistically isotropic noise)

  Args:
    :n (*int*): Number of maps
    :npix (*int*): Number of pixel
    :lmax (*int*): Maximum multipole of alm
    :clh[*n,l*] (*double*): Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
    :nij[*n,pix*] (*double*): Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
    :b[*n,l,m*] (*dcmplx*): RHS, with bounds (0:n-1,0:lmax,0:lmax)
    :itern (*int*): Number of interation
    
  Args(optional): 
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6

  Returns:
    :x[*n,l,m*] (*dcmplx*): C-inverse filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
 

  Usage:
    :x = curvedsky.cninv.cg_algorithm(n,npix,lmax,clh,nij,b,itern,eps):
  """
  return libcurvedsky.cninv.cg_algorithm(n,npix,lmax,clh,nij,b,itern,eps)

def mat_multi(n,npix,lmax,clh,nij,x,mtype):
  """
 multiplying matrix
  Usage:
    :v = curvedsky.cninv.mat_multi(n,npix,lmax,clh,nij,x,mtype):
  """
  return libcurvedsky.cninv.mat_multi(n,npix,lmax,clh,nij,x,mtype)

