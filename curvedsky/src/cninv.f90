!////////////////////////////////////////////////////!
! Map filtering
!////////////////////////////////////////////////////!

module cninv
  !from F90/src_utils
  use constants, only: pi
  use spht
  implicit none

  private pi

contains


subroutine cnfilter(npix,lmax,cl,nij,alm,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*    1) The signal power spectrum is isotropic Gaussian. 
!*    2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*    3) The data model is bxS+N
!*
!*  Args:
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:lmax)
!*    :nij[pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:npix-1)
!*    :alm[l,m] (dcmplx) : Input alm, with bouds (0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*
!*  Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!*  Returns:
!*    :xlm[l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  character(1), intent(in) :: filter
  integer, intent(in) :: npix, lmax, itern
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  !opt4py :: filter = ''
  double precision, intent(in), dimension(0:lmax) :: cl
  double precision, intent(in), dimension(0:npix-1) :: nij
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: xlm
  !internal
  integer :: ni, l
  double precision :: clh(1:1,0:lmax,0:lmax), ninv(1:1,0:npix-1)
  double complex :: b(1:1,0:lmax,0:lmax), nalm(1:1,0:lmax,0:lmax), nxlm(1:1,0:lmax,0:lmax)

  clh = 0d0
  do l = 1, lmax
    clh(1,l,0:l) = dsqrt(cl(l))
  end do

  !first compute b = C^1/2 N^-1 X
  ninv(1,:)   = nij
  nalm(1,:,:) = alm
  call mat_multi(1,npix,lmax,clh,ninv,nalm,b,'rhs')
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm(1,npix,lmax,clh,ninv,b,nxlm,itern,eps)

  xlm = 0d0
  do l = 1, lmax
      if (filter=='')   xlm(l,0:l) = nxlm(1,l,0:l)/clh(1,l,0:l)
      if (filter=='W')  xlm(l,0:l) = nxlm(1,l,0:l)*clh(1,l,0:l)
  end do

end subroutine cnfilter


subroutine cnfilterpol(n,npix,lmax,cl,nij,alm,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is bxS+N
!*
!*  Args:
!*    :n (int) : Number of maps
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :nij[n,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :alm[n,l,m] (dcmplx) : Input alm, with bouds (0:n-1,0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*
!*  Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!*  Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  character(1), intent(in) :: filter
  integer, intent(in) :: n, npix, lmax, itern
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  !opt4py :: filter = ''
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  integer :: ni, l
  double precision, dimension(n,0:lmax,0:lmax) :: clh
  double complex :: b(n,0:lmax,0:lmax)

  clh = 0d0
  do ni = 1, n
    do l = 1, lmax
      clh(ni,l,0:l) = dsqrt(cl(ni,l))
    end do
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi(n,npix,lmax,clh,nij,alm,b,'rhs')
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm(n,npix,lmax,clh,nij,b,xlm,itern,eps)

  do l = 1, lmax
    if (filter=='')   xlm(:,l,0:l) = xlm(:,l,0:l)/clh(:,l,0:l)
    if (filter=='W')  xlm(:,l,0:l) = xlm(:,l,0:l)*clh(:,l,0:l)
  end do

end subroutine cnfilterpol


subroutine cg_algorithm(n,npix,lmax,clh,nij,b,x,itern,eps)
!*  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
!*  The code assumes
!*    1) A = [1 + C^1/2 N^-1 C^1/2]
!*    2) C^1/2 is diagonal
!*    3) N is diagonal in pixel space (statistically isotropic noise)
!*
!*  Args:
!*    :n (int) : Number of maps
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :clh[n,l] (double) : Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
!*    :nij[n,pix] (double) : Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
!*    :b[n,l,m] (dcmplx) : RHS, with bounds (0:n-1,0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*    
!*  Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*
!*  Returns:
!*    :x[n,l,m] (dcmplx) : C-inverse filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!* 
!*
  implicit none
  !I/O
  integer, intent(in) :: n, npix, lmax, itern
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  double precision, intent(in), dimension(n,0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: b
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: x
  !internal
  integer :: ni, i, l, ro=50
  double precision :: absb, absr, d, d0, td, alpha, ndiag, M(1:n,0:lmax,0:lmax)
  double complex, dimension(1:n,0:lmax,0:lmax) :: r, z, p, Ap

  !preconditioning matrix (M) which makes MA ~ I
  do ni = 1, n
    ndiag = sum(nij(ni,:))*4d0*pi/size(nij)
    M(ni,:,:) = 1d0/(1d0+ndiag*clh(ni,:,:)**2)
  end do
  absb = dsqrt(sum(abs(b)**2))

  !initial value (this is the solution if MA=I)
  x = M*b

  !residual
  call mat_multi(n,npix,lmax,clh,nij,x,r,'lhs')
  r = b - r

  !set other values
  z = M*r 
  p = z

  !initial distance
  d0 = sum(conjg(r)*z)
  d  = d0

  do i = 1, itern

    call mat_multi(n,npix,lmax,clh,nij,p,Ap,'lhs')
    alpha = d/sum(conjg(p)*Ap)
    x = x + alpha*p
    r = r - alpha*Ap

    absr = dsqrt(sum(abs(r)**2))
    if (i-int(dble(i)/dble(ro))*ro==0)  write(*,*) absr/absb, d/d0

    z = M*r
    td = sum(conjg(r)*z)
    p  = z + (td/d)*p
    d  = td

    ! check exit condition
    if (absr<eps*absb) then 
      !exit loop if |r|/|b| becomes very small
      write(*,*) i, absr/absb
      exit !Norm of r is sufficiently small
    end if

  end do


end subroutine cg_algorithm



subroutine cnfilter_lat(n,k1,lmax,cl,bl1,npix1,nij1,map1,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is bxS+N
!*
!*  Args:
!*    :n (int) : T(1), Q/U(2) or T/Q/U(3)
!*    :k1 (int) : Number of frequencies
!*    :npix1 (int) : Number of pixels for each input maps and inv noise covariance
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :bl1[k,l] (double) : Beam spectrum, with bounds (0:k1-1,0:lmax)
!*    :nij1[n,k,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k1-1,0:npix1-1)
!*    :map1[n,k,pix] (double) : Input maps, with bouds (0:n-1,0:k1-1,0:npix1-1)
!*    :itern (int) : Number of interation
!*
!*  Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!*  Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  character(1), intent(in) :: filter
  integer, intent(in) :: n, k1, npix1, lmax, itern
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  !opt4py :: filter = ''
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(k1,0:lmax) :: bl1
  double precision, intent(in), dimension(n,k1,0:npix1-1) :: nij1, map1
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  integer :: ni, l, ki
  double precision, dimension(n,k1,0:lmax,0:lmax) :: clh1
  double complex, dimension(n,0:lmax,0:lmax) :: b1

  clh1 = 0d0
  do ni = 1, n
    do l = 1, lmax
      do ki = 1, k1
        clh1(ni,ki,l,0:l) = dsqrt(cl(ni,l))*bl1(ki,l)
      end do
    end do
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi_freq_rhs(n,k1,npix1,lmax,clh1,nij1,map1,b1)
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm_lat(n,k1,lmax,clh1,npix1,nij1,b1,xlm,itern,eps)

  do l = 1, lmax
    do ni = 1, n
      if (filter=='')   xlm(ni,l,0:l) = xlm(ni,l,0:l)/dsqrt(cl(ni,l))
      if (filter=='W')  xlm(ni,l,0:l) = xlm(ni,l,0:l)*dsqrt(cl(ni,l))
    end do
  end do

end subroutine cnfilter_lat


subroutine cg_algorithm_lat(n,k1,lmax,clh1,npix1,nij1,b,x,itern,eps)
!*  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
!*  The code assumes 
!*    1) A = [1 + C^1/2 N^-1 C^1/2]
!*    2) C^1/2 is diagonal
!*    3) N is diagonal in pixel space (statistically isotropic noise)
!*
!*  Args:
!*    :n (int) : T(1), Q/U(2), or T/Q/U(3)
!*    :k1 (int) : Number of freq
!*    :npix1 (int) : Number of pixels
!*    :lmax (int) : Maximum multipole of alm
!*    :clh1[n,k,l] (double) : Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
!*    :nij1[n,k,pix] (double) : Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
!*    :b[n,l,m] (dcmplx) : RHS, with bounds (0:n-1,0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*    
!*  Args(optional):
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*
!*  Returns:
!*    :x[n,l,m] (dcmplx) : C-inverse filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!* 
!*
  implicit none
  !I/O
  integer, intent(in) :: n, k1, lmax, itern, npix1
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  double precision, intent(in), dimension(n,k1,0:lmax,0:lmax) :: clh1
  double precision, intent(in), dimension(n,k1,0:npix1-1) :: nij1
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: b
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: x
  !internal
  integer :: ni, ki, i, l, ro=50
  double precision :: absb, absr, d, d0, td, alpha, M(1:n,0:lmax,0:lmax), mm(0:lmax,0:lmax)
  double complex, dimension(1:n,0:lmax,0:lmax) :: r, z, p, Ap

  !preconditioning matrix (M) which makes MA ~ I
  do ni = 1, n
    mm = 0d0
    do ki = 1, k1
      mm = mm + clh1(ni,ki,:,:)**2 * sum(nij1(ni,ki,:))*4d0*pi/size(nij1(:,ki,:))
    end do
    M(ni,:,:) = 1d0/(1d0+mm)
  end do
  absb = dsqrt(sum(abs(b)**2))

  !initial value (this is the solution if MA=I)
  x = M*b

  !residual
  call mat_multi_freq_lhs(n,k1,npix1,lmax,clh1,nij1,x,r)
  r = b - r

  !set other values
  z = M*r
  p = z

  !initial distance
  d0 = sum(conjg(r)*z)
  d  = d0

  do i = 1, itern

    call mat_multi_freq_lhs(n,k1,npix1,lmax,clh1,nij1,p,Ap)
    alpha = d/sum(conjg(p)*Ap)
    x = x + alpha*p
    r = r - alpha*Ap

    absr = dsqrt(sum(abs(r)**2))
    if (i-int(dble(i)/dble(ro))*ro==0)  write(*,*) absr/absb, d/d0

    z = M*r
    td = sum(conjg(r)*z)
    p  = z + (td/d)*p
    d  = td

    ! check exit condition
    if (absr<eps*absb) then 
      !exit loop if |r|/|b| becomes very small
      write(*,*) i, absr/absb
      exit !Norm of r is sufficiently small
    end if

  end do

end subroutine cg_algorithm_lat


subroutine cnfilter_so(n,k1,k2,lmax,cl,bl1,bl2,npix1,npix2,nij1,nij2,map1,map2,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is bxS+N
!*
!*  Args:
!*    :n (int) : T(1), Q/U(2) or T/Q/U(3)
!*    :k1 (int) : Number of frequencies
!*    :k2 (int) : Number of frequencies
!*    :npix1 (int) : Number of pixels for each input maps and inv noise covariance
!*    :npix2 (int) : Number of pixels for each input maps and inv noise covariance
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :bl1[k,l] (double) : Beam spectrum, with bounds (0:k1-1,0:lmax)
!*    :bl2[k,l] (double) : Beam spectrum, with bounds (0:k2-1,0:lmax)
!*    :nij1[n,k,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k1-1,0:npix1-1)
!*    :nij2[n,k,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k2-1,0:npix2-1)
!*    :map1[n,k,pix] (double) : Input maps, with bouds (0:n-1,0:k1-1,0:npix1-1)
!*    :map2[n,k,pix] (double) : Input maps, with bouds (0:n-1,0:k2-1,0:npix2-1)
!*    :itern (int) : Number of interation
!*
!*  Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!*  Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  character(1), intent(in) :: filter
  integer, intent(in) :: n, k1, k2, npix1, npix2, lmax, itern
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  !opt4py :: filter = ''
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(k1,0:lmax) :: bl1
  double precision, intent(in), dimension(k2,0:lmax) :: bl2
  double precision, intent(in), dimension(n,k1,0:npix1-1) :: nij1, map1
  double precision, intent(in), dimension(n,k2,0:npix2-1) :: nij2, map2
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  integer :: ni, l, ki
  double precision, dimension(n,k1,0:lmax,0:lmax) :: clh1
  double precision, dimension(n,k2,0:lmax,0:lmax) :: clh2
  double complex, dimension(n,0:lmax,0:lmax) :: b1, b2

  clh1 = 0d0
  clh2 = 0d0
  do ni = 1, n
    do l = 1, lmax
      do ki = 1, k1
        clh1(ni,ki,l,0:l) = dsqrt(cl(ni,l))*bl1(ki,l)
      end do
      do ki = 1, k2
        clh2(ni,ki,l,0:l) = dsqrt(cl(ni,l))*bl2(ki,l)
      end do
    end do
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi_freq_rhs(n,k1,npix1,lmax,clh1,nij1,map1,b1)
  call mat_multi_freq_rhs(n,k2,npix2,lmax,clh2,nij2,map2,b2)
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm_so(n,k1,k2,lmax,clh1,clh2,npix1,npix2,nij1,nij2,b1+b2,xlm,itern,eps)

  do l = 1, lmax
    do ni = 1, n
      if (filter=='')   xlm(ni,l,0:l) = xlm(ni,l,0:l)/dsqrt(cl(ni,l))
      if (filter=='W')  xlm(ni,l,0:l) = xlm(ni,l,0:l)*dsqrt(cl(ni,l))
    end do
  end do

end subroutine cnfilter_so


subroutine cg_algorithm_so(n,k1,k2,lmax,clh1,clh2,npix1,npix2,nij1,nij2,b,x,itern,eps)
!*  Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
!*  The code assumes 
!*    1) A = [1 + C^1/2 N^-1 C^1/2]
!*    2) C^1/2 is diagonal
!*    3) N is diagonal in pixel space (statistically isotropic noise)
!*
!*  Args:
!*    :n (int) : T(1), Q/U(2), or T/Q/U(3)
!*    :k1 (int) : Number of freq
!*    :k2 (int) : Number of freq
!*    :npix1 (int) : Number of pixels
!*    :npix2 (int) : Number of pixels
!*    :lmax (int) : Maximum multipole of alm
!*    :clh1[n,k,l] (double) : Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
!*    :clh2[n,k,l] (double) : Square root of angular spectrum (C^1/2), with bounds (0:n-1,0:lmax)
!*    :nij1[n,k,pix] (double) : Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
!*    :nij2[n,k,pix] (double) : Inverse of the noise variance (N^-1) at each pixel, with bounds (0:n-1,0:npix-1)
!*    :b[n,l,m] (dcmplx) : RHS, with bounds (0:n-1,0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*    
!*  Args(optional):
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*
!*  Returns:
!*    :x[n,l,m] (dcmplx) : C-inverse filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!* 
!*
  implicit none
  !I/O
  integer, intent(in) :: n, k1, k2, lmax, itern, npix1, npix2
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  double precision, intent(in), dimension(n,k1,0:lmax,0:lmax) :: clh1
  double precision, intent(in), dimension(n,k2,0:lmax,0:lmax) :: clh2
  double precision, intent(in), dimension(n,k1,0:npix1-1) :: nij1
  double precision, intent(in), dimension(n,k2,0:npix2-1) :: nij2
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: b
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: x
  !internal
  integer :: ni, ki, i, l, ro=50
  double precision :: absb, absr, d, d0, td, alpha, M(1:n,0:lmax,0:lmax), mm(0:lmax,0:lmax)
  double complex, dimension(1:n,0:lmax,0:lmax) :: r, r1, r2, z, p, Ap, Ap1, Ap2

  !preconditioning matrix (M) which makes MA ~ I
  do ni = 1, n
    mm = 0d0
    do ki = 1, k1
      mm = mm + clh1(ni,ki,:,:)**2 * sum(nij1(ni,ki,:))*4d0*pi/size(nij1(:,ki,:))
    end do
    do ki = 1, k2
      mm = mm + clh2(ni,ki,:,:)**2 * sum(nij2(ni,ki,:))*4d0*pi/size(nij2(:,ki,:))
    end do
    M(ni,:,:) = 1d0/(1d0+mm)
  end do
  absb = dsqrt(sum(abs(b)**2))

  !initial value (this is the solution if MA=I)
  x = M*b

  !residual
  call mat_multi_freq_lhs(n,k1,npix1,lmax,clh1,nij1,x,r1)
  call mat_multi_freq_lhs(n,k2,npix2,lmax,clh2,nij2,x,r2)
  r = r1+r2-x
  r = b - r

  !set other values
  z = M*r
  p = z

  !initial distance
  d0 = sum(conjg(r)*z)
  d  = d0

  do i = 1, itern

    call mat_multi_freq_lhs(n,k1,npix1,lmax,clh1,nij1,p,Ap1)
    call mat_multi_freq_lhs(n,k2,npix2,lmax,clh2,nij2,p,Ap2)
    Ap = Ap1+Ap2-p
    alpha = d/sum(conjg(p)*Ap)
    x = x + alpha*p
    r = r - alpha*Ap

    absr = dsqrt(sum(abs(r)**2))
    if (i-int(dble(i)/dble(ro))*ro==0)  write(*,*) absr/absb, d/d0

    z = M*r
    td = sum(conjg(r)*z)
    p  = z + (td/d)*p
    d  = td

    ! check exit condition
    if (absr<eps*absb) then 
      !exit loop if |r|/|b| becomes very small
      write(*,*) i, absr/absb
      exit !Norm of r is sufficiently small
    end if

  end do

end subroutine cg_algorithm_so


end module cninv


