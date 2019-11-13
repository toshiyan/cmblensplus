import libcurvedsky

def cnfilter(npix,lmax,cl,nij,alm,itern,eps=1e-6,filter=''):
  """
 Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
 This code assumes
    1) The signal power spectrum is isotropic Gaussian. 
    2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
    3) The data model is bxS+N

  Args:
    :npix (*int*): Number of pixel
    :lmax (*int*): Maximum multipole of alm
    :cl[*n,l*] (*double*): Angular power spectrum of alm, with bounds (0:lmax)
    :nij[*pix*] (*double*): Inverse of the noise variance at each pixel, with bounds (0:npix-1)
    :alm[*l,m*] (*dcmplx*): Input alm, with bouds (0:lmax,0:lmax)
    :itern (*int*): Number of interation

  Args(optional): 
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
    :filter (*str*): C-inverse ('') or Wiener filter (W), default to C-inverse.

  Returns:
    :xlm[*l,m*] (*dcmplx*): C-inverse / Wiener filtered multipoles, with bounds (0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.cninv.cnfilter(npix,lmax,cl,nij,alm,itern,eps,filter):
  """
  return libcurvedsky.cninv.cnfilter(npix,lmax,cl,nij,alm,itern,eps,filter)

def cnfilterpol(n,npix,lmax,cl,nij,alm,itern,eps=1e-6,filter=''):
  """
 Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
 This code assumes
   1) The signal power spectrum is isotropic Gaussian. 
   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
   3) The data model is bxS+N

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
    :xlm = curvedsky.cninv.cnfilterpol(n,npix,lmax,cl,nij,alm,itern,eps,filter):
  """
  return libcurvedsky.cninv.cnfilterpol(n,npix,lmax,cl,nij,alm,itern,eps,filter)

def cg_algorithm(n,npix,lmax,clh,nij,b,itern,eps=1e-6):
  """
  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
  The code assumes
    1) A = [1 + C^1/2 N^-1 C^1/2]
    2) C^1/2 is diagonal
    3) N is diagonal in pixel space (statistically isotropic noise)

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

def cnfilter_lat(n,k1,lmax,cl,bl1,npix1,nij1,map1,itern,eps=1e-6,filter=''):
  """
 Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
 This code assumes
   1) The signal power spectrum is isotropic Gaussian. 
   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
   3) The data model is bxS+N

  Args:
    :n (*int*): T(1), Q/U(2) or T/Q/U(3)
    :k1 (*int*): Number of frequencies
    :npix1 (*int*): Number of pixels for each input maps and inv noise covariance
    :lmax (*int*): Maximum multipole of alm
    :cl[*n,l*] (*double*): Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
    :bl1[*k,l*] (*double*): Beam spectrum, with bounds (0:k1-1,0:lmax)
    :nij1[*n,k,pix*] (*double*): Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k1-1,0:npix1-1)
    :map1[*n,k,pix*] (*double*): Input maps, with bouds (0:n-1,0:k1-1,0:npix1-1)
    :itern (*int*): Number of interation

  Args(optional): 
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
    :filter (*str*): C-inverse ('') or Wiener filter (W), default to C-inverse.

  Returns:
    :xlm[*n,l,m*] (*dcmplx*): C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.cninv.cnfilter_lat(n,k1,lmax,cl,bl1,npix1,nij1,map1,itern,eps,filter):
  """
  return libcurvedsky.cninv.cnfilter_lat(n,k1,lmax,cl,bl1,npix1,nij1,map1,itern,eps,filter)

def cg_algorithm_lat(n,k1,lmax,clh1,npix1,nij1,b,itern,eps=1e-6):
  """
  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
  The code assumes 
    1) A = [1 + C^1/2 N^-1 C^1/2]
    2) C^1/2 is diagonal
    3) N is diagonal in pixel space (statistically isotropic noise)

  Args:
    :n (*int*): T(1), Q/U(2), or T/Q/U(3)
    :k1 (*int*): Number of freq
    :npix1 (*int*): Number of pixels
    :lmax (*int*): Maximum multipole of alm
    :clh1[*n,k,l*] (*double*): Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
    :nij1[*n,k,pix*] (*double*): Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
    :b[*n,l,m*] (*dcmplx*): RHS, with bounds (0:n-1,0:lmax,0:lmax)
    :itern (*int*): Number of interation
    
  Args(optional):
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6

  Returns:
    :x[*n,l,m*] (*dcmplx*): C-inverse filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
 

  Usage:
    :x = curvedsky.cninv.cg_algorithm_lat(n,k1,lmax,clh1,npix1,nij1,b,itern,eps):
  """
  return libcurvedsky.cninv.cg_algorithm_lat(n,k1,lmax,clh1,npix1,nij1,b,itern,eps)

def cnfilter_so(n,k1,k2,lmax,cl,bl1,bl2,npix1,npix2,nij1,nij2,map1,map2,itern,eps=1e-6,filter=''):
  """
 Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
 This code assumes
   1) The signal power spectrum is isotropic Gaussian. 
   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
   3) The data model is bxS+N

  Args:
    :n (*int*): T(1), Q/U(2) or T/Q/U(3)
    :k1 (*int*): Number of frequencies
    :k2 (*int*): Number of frequencies
    :npix1 (*int*): Number of pixels for each input maps and inv noise covariance
    :npix2 (*int*): Number of pixels for each input maps and inv noise covariance
    :lmax (*int*): Maximum multipole of alm
    :cl[*n,l*] (*double*): Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
    :bl1[*k,l*] (*double*): Beam spectrum, with bounds (0:k1-1,0:lmax)
    :bl2[*k,l*] (*double*): Beam spectrum, with bounds (0:k2-1,0:lmax)
    :nij1[*n,k,pix*] (*double*): Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k1-1,0:npix1-1)
    :nij2[*n,k,pix*] (*double*): Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k2-1,0:npix2-1)
    :map1[*n,k,pix*] (*double*): Input maps, with bouds (0:n-1,0:k1-1,0:npix1-1)
    :map2[*n,k,pix*] (*double*): Input maps, with bouds (0:n-1,0:k2-1,0:npix2-1)
    :itern (*int*): Number of interation

  Args(optional): 
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
    :filter (*str*): C-inverse ('') or Wiener filter (W), default to C-inverse.

  Returns:
    :xlm[*n,l,m*] (*dcmplx*): C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)

  Usage:
    :xlm = curvedsky.cninv.cnfilter_so(n,k1,k2,lmax,cl,bl1,bl2,npix1,npix2,nij1,nij2,map1,map2,itern,eps,filter):
  """
  return libcurvedsky.cninv.cnfilter_so(n,k1,k2,lmax,cl,bl1,bl2,npix1,npix2,nij1,nij2,map1,map2,itern,eps,filter)

def cg_algorithm_so(n,k1,k2,lmax,clh1,clh2,npix1,npix2,nij1,nij2,b,itern,eps=1e-6):
  """
  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
  The code assumes 
    1) A = [1 + C^1/2 N^-1 C^1/2]
    2) C^1/2 is diagonal
    3) N is diagonal in pixel space (statistically isotropic noise)

  Args:
    :n (*int*): T(1), Q/U(2), or T/Q/U(3)
    :k1 (*int*): Number of freq
    :k2 (*int*): Number of freq
    :npix1 (*int*): Number of pixels
    :npix2 (*int*): Number of pixels
    :lmax (*int*): Maximum multipole of alm
    :clh1[*n,k,l*] (*double*): Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
    :clh2[*n,k,l*] (*double*): Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
    :nij1[*n,k,pix*] (*double*): Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
    :nij2[*n,k,pix*] (*double*): Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
    :b[*n,l,m*] (*dcmplx*): RHS, with bounds (0:n-1,0:lmax,0:lmax)
    :itern (*int*): Number of interation
    
  Args(optional):
    :eps (*double*): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6

  Returns:
    :x[*n,l,m*] (*dcmplx*): C-inverse filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
 

  Usage:
    :x = curvedsky.cninv.cg_algorithm_so(n,k1,k2,lmax,clh1,clh2,npix1,npix2,nij1,nij2,b,itern,eps):
  """
  return libcurvedsky.cninv.cg_algorithm_so(n,k1,k2,lmax,clh1,clh2,npix1,npix2,nij1,nij2,b,itern,eps)

