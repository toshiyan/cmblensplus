!////////////////////////////////////////////////////!
! Map filtering
!////////////////////////////////////////////////////!

module cninv
  !from F90/src_utils
  use spht,      only: spht_alm2map, spht_map2alm
  implicit none

  private spht_alm2map, spht_map2alm

contains


subroutine cnfilter(n,npix,lmax,cl,nij,alm,itern,xlm,eps,filter)
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
  double precision, dimension(n,0:lmax) :: clh
  double complex :: b(n,0:lmax,0:lmax)

  clh = 0d0
  do l = 1, lmax
    clh(:,l) = dsqrt(cl(:,l))
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi(n,npix,lmax,clh,nij,alm,b,'rhs')
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm(n,npix,lmax,clh,nij,b,xlm,itern,eps)

  do ni = 1, n
    do l = 1, lmax
      if (filter=='')   xlm(ni,l,:) = xlm(ni,l,:)/clh(ni,l)
      if (filter=='W')  xlm(ni,l,:) = xlm(ni,l,:)*clh(ni,l)
    end do
  end do

end subroutine cnfilter


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
  double precision, intent(in), dimension(n,0:lmax) :: clh
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: b
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: x
  !internal
  integer :: ni, i, l, iter, roundoff=50
  double precision :: dr, dr0, tdr, alpha, ndiag
  double complex, dimension(1:n,0:lmax,0:lmax) :: r, z, p, r0, Ap, M

  !preconditioning matrix (M) which makes MA ~ I
  M = 0d0
  do ni = 1, n
    ndiag = sum(nij(ni,:))/size(nij)
    do l = 1, lmax
      M(ni,l,:) = 1d0/(1d0+ndiag*clh(ni,l)**2)
    end do
  end do

  !initial value (this is the solution if MA=I)
  x = M*b
  !x = 0

  !residual
  call mat_multi(n,npix,lmax,clh,nij,x,r,'lhs')
  r = b - r 

  !set other values
  z = M*r 
  p = z

  !initial distance
  dr0 = sum(conjg(r)*z)
  dr  = dr0
  if (dr<1d-10*size(x)) return 

  !iteration num
  iter = 0

  do i = 1, itern

    ! check exit condition
    if (dr<eps*size(x)) then 
      !exit loop if |r_i|/|r_0| becomes very small
      write(*,*) iter, dr/dr0
      exit !Norm of r is sufficiently small
    end if
    iter = iter + 1

    ! Ap(i)
    call mat_multi(n,npix,lmax,clh,nij,p,Ap,'lhs')
    
    ! alpha(i) = |r(i)|^2_tA / |p(i)|^2_A
    alpha = dr/sum(conjg(p)*Ap)

    ! x(i+1) = x(i) + akpha(i)p(i)
    x = x + alpha*p

    ! r(i+1) = r(i) - alpha(i)Ap(i)
    r = r - alpha*Ap

    ! residual multiplied by preconditioning matrix
    z = M*r

    !distance
    tdr = sum(conjg(r)*z)

    ! p(i+1) = z(i+1) + beta(i)p(i)
    p   = z + (tdr/dr)*p

    ! update dr
    dr  = tdr

    !check distance
    if(iter-int(dble(iter)/dble(roundoff))*roundoff==0) then 
      !r0 = r
      !call mat_multi(npix,lmax,clh,nij,x,r0,'lhs')  !r0 = A x
      !r0 = b - r0
      write(*,*) iter, dr/dr0
    end if

  end do


end subroutine cg_algorithm


subroutine mat_multi(n,npix,lmax,clh,nij,x,v,mtype)
!* multiplying matrix
  implicit none
  !I/O
  character(3), intent(in) :: mtype
  integer, intent(in) :: n, npix, lmax
  double precision, intent(in), dimension(n,0:lmax) :: clh
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: x
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  integer :: p, l, nside
  double precision :: map(n,0:npix-1)
  double complex :: alm(n,0:lmax,0:lmax)

  nside = int(sqrt(npix/12d0))

  if (mtype=='lhs') then
    do p = 1, n
      do l = 0, lmax
        !multiply C^{1/2}
        alm(p,l,:) = clh(p,l)*x(p,l,:)
      end do
    end do
  else
    alm = x
  end if

  !multiply noise covariance in pixel space
  call spht_alm2map(n,npix,lmax,lmax,alm,map)
  map = map*nij
  call spht_map2alm(n,npix,lmax,lmax,map,alm)

  !multiply C^{1/2}
  do p = 1, n
    do l = 0, lmax
      alm(p,l,:) = clh(p,l)*alm(p,l,:)
    end do
  end do

  if (mtype=='lhs') then
    !make 1+C^{1/2}N^{-1}C^{1/2}
    v = x + alm
  else
    v = alm
  end if

end subroutine mat_multi


end module cninv



