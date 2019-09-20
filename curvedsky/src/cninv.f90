!////////////////////////////////////////////////////!
! Map filtering
!////////////////////////////////////////////////////!

module cninv
  !from F90/src_utils
  use constants, only: pi
  use spht,      only: mat_multi, mat_multi0
  implicit none

  private pi
  private mat_multi, mat_multi0

contains


subroutine cnfilter(npix,lmax,cl,nij,alm,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   - The signal power spectrum is isotropic Gaussian. 
!*   - Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   - The data model is bxS+N
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


subroutine cnfilter0(npix,lmax,cl,nij,alm,itern,xlm)
!* Computing inverse-variance filtered multipoles: C^-1d
!* This code assumes
!*   - The signal power spectrum is isotropic Gaussian. 
!*   - Inverse noise covariance is given in pixel space and uncorrelated (nij = sigma x delta_ij).
!*   - The data model is bxS+N
!*
!*  Args:
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[l] (double) : Angular power spectrum of alm, with bounds (0:lmax)
!*    :nij[pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:npix-1)
!*    :alm[l,m] (dcmplx) : Input alm, with bouds (0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*
!*  Returns:
!*    :xlm[l,m] (dcmplx) : C-inverse filtered multipoles, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, itern
  double precision, intent(in), dimension(0:lmax) :: cl
  double precision, intent(in), dimension(0:npix-1) :: nij
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: xlm
  !internal
  integer :: l
  double precision, dimension(0:lmax,0:lmax) :: clh
  double complex, allocatable :: b(:,:)

  clh = 0d0
  do l = 1, lmax
    clh(l,0:l) = dsqrt(cl(l))
  end do

  allocate(b(0:lmax,0:lmax));  b=0d0
  !first compute b = C^1/2 N^-1 X
  call mat_multi0(npix,lmax,clh,nij,alm,b)
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm0(npix,lmax,clh,nij,b,xlm,itern)
  deallocate(b)

  do l = 1, lmax
    xlm(l,0:l) = xlm(l,0:l)*clh(l,0:l)
  end do

end subroutine cnfilter0


subroutine cnfilterpol(n,npix,lmax,cl,nij,alm,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   - The signal power spectrum is isotropic Gaussian. 
!*   - Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   - The data model is bxS+N
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
!*    - A = [1 + C^1/2 N^-1 C^1/2]
!*    - C^1/2 is diagonal
!*    - N is diagonal in pixel space (statistically isotropic noise)
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
  double precision :: bb, d, d0, td, alpha, ndiag, M(1:n,0:lmax,0:lmax)
  double complex, dimension(1:n,0:lmax,0:lmax) :: r, z, p, Ap

  !preconditioning matrix (M) which makes MA ~ I
  do ni = 1, n
    ndiag = sum(nij(ni,:))*4d0*pi/size(nij)
    M(ni,:,:) = 1d0/(1d0+ndiag*clh(ni,:,:)**2)
  end do
  !write(*,*) sum(M)/size(M)
  bb = sum(abs(b)**2)

  !initial value (this is the solution if MA=I)
  x = M*b
  !write(*,*) sum(abs(x)**2)/bb, bb
  !x = 0

  !residual
  call mat_multi(n,npix,lmax,clh,nij,x,r,'lhs')
  !write(*,*) sum(abs(r)**2)/bb
  r = b - r
  !write(*,*) sum(abs(r)**2)/bb

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

    if (i-int(dble(i)/dble(ro))*ro==0)  write(*,*) sum(abs(r)**2)/bb, d/bb

    z = M*r
    td = sum(conjg(r)*z)
    p  = z + (td/d)*p
    d  = td

    ! check exit condition
    if (d<eps**2*d0) then 
      !exit loop if |r_i|/|r_0| becomes very small
      write(*,*) i, d/d0
      exit !Norm of r is sufficiently small
    end if

  end do


end subroutine cg_algorithm


subroutine cg_algorithm0(npix,lmax,clh,nij,b,x,itern)
!* Searching a solution of Ax=b with the Conjugate Gradient iteratively
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax, itern
  double precision, intent(in), dimension(0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(0:npix-1) :: nij
  double complex, intent(in), dimension(0:lmax,0:lmax) :: b
  double complex, intent(out), dimension(0:lmax,0:lmax) :: x
  !internal
  integer :: i, l, ro=50
  double precision :: bb, d, d0, td, alpha, eps = 1d-6, ndiag, M(0:lmax,0:lmax)
  double complex, dimension(0:lmax,0:lmax) :: r, tAr, p, r0, Ap, z

  ndiag = sum(nij)*4d0*pi/dble(npix) !diagonal inverse noise covariance
  M = 1d0/(1d0+ndiag*clh**2) !preconditioning matrix
  !write(*,*) sum(M)/size(M)
  bb = sum(abs(b)**2)

  x = M*b
  !write(*,*) sum(abs(x)**2)/bb, bb
  !x = 0

  call mat_multi0(npix,lmax,clh,nij,x,r,'lhs')  ! r = Ax
  !write(*,*) sum(abs(r)**2)/bb
  r = b - r
  !write(*,*) sum(abs(r)**2)/bb

  !multiply preconditioner
  z = M*r
  p = z

  !initial distance
  d0 = sum(conjg(r)*z)
  d  = d0

  do i = 1, itern

    call mat_multi0(npix,lmax,clh,nij,p,Ap,'lhs')  ! Ap = A p
    alpha = d/sum(conjg(p)*Ap)
    x = x + alpha*p
    r = r - alpha*Ap

    if (i-int(dble(i)/dble(ro))*ro==0)  write(*,*) sum(abs(r)**2)/bb, d/bb

    z  = M*r
    td = sum(conjg(r)*z)
    p  = z + (td/d)*p
    d  = td

    if (d<eps**2*d0) then 
      write(*,*) i, d/d0
      exit !Norm of r is sufficiently small
    end if 

  end do


end subroutine cg_algorithm0



end module cninv



