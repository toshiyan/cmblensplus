!////////////////////////////////////////////////////!
! Conjugate Gradient Decent
!////////////////////////////////////////////////////!

module hp_cgd
  ! from src_utils
  use constants, only: pi, iu
  use general, only: ave, str
  ! from src_matrix
  use lapack95, only: inv_matrix_c
  ! from modules in this directory
  use hp_udgrade, only: udgrade_ring_1d_d
  use hp_spht, only: spht_alm2map, spht_map2alm
  use hp_utils, only: trans_alm2array_1d, trans_array2alm_1d, trans_array2alm_2d
  implicit none

  ! for multi-grid method
  type mg_covmat
    double precision, allocatable :: nij(:,:), clh(:,:,:)
  end type mg_covmat

  type mg_chain
    type(mg_covmat), allocatable :: cv(:,:)
    logical :: verbose
    integer :: n, lsp, matn, ro
    integer, allocatable :: mnmax(:), npix(:,:), nside(:,:), lmax(:), itn(:)
    double precision, allocatable :: eps(:), pmat(:,:,:)
    double complex, allocatable :: mat(:,:), imat(:,:)
  end type mg_chain

  private pi, iu
  private spht_alm2map, spht_map2alm
  private inv_matrix_c
  private udgrade_ring_1d_d

contains


subroutine set_mgchain(mgc,chn,mn,mnmaxs,lmaxs,nsides,itns,eps,verbose,ro)
  implicit none
  type(mg_chain), intent(out) :: mgc
  logical :: verbose
  integer, intent(in) :: chn, mn, ro, mnmaxs(chn), itns(chn), lmaxs(chn), nsides(chn,mn)
  double precision, intent(in) :: eps(chn)

  mgc%n = chn
  allocate(mgc%mnmax(chn),mgc%nside(chn,mn),mgc%lmax(chn),mgc%itn(chn),mgc%eps(chn),mgc%npix(chn,mn))

  mgc%nside = nsides
  mgc%lmax  = lmaxs
  mgc%itn   = itns
  mgc%eps   = eps
  mgc%npix  = 12*nsides**2
  mgc%verbose = verbose
  mgc%lsp   = lmaxs(chn)
  mgc%mnmax = mnmaxs
  mgc%ro    = ro

  write(*,'('//str(chn)//'(I4,X))') mgc%nside
  write(*,'('//str(chn)//'(I4,X))') mgc%mnmax
  write(*,'('//str(chn)//'(I4,X))') mgc%lmax
  write(*,'('//str(chn)//'(I4,X))') mgc%itn
  write(*,'('//str(chn)//'(E14.7,X))') mgc%eps
  write(*,*) mgc%lsp
  if (mgc%verbose) write(*,*) mgc%verbose
  if (mgc%ro/=50)  write(*,*) mgc%ro

end subroutine set_mgchain


subroutine free_mgchain(mgc)
  implicit none
  type(mg_chain), intent(inout) :: mgc

  deallocate(mgc%npix,mgc%nside,mgc%itn,mgc%lmax,mgc%eps)
  if (allocated(mgc%imat)) deallocate(mgc%imat)
  if (allocated(mgc%pmat)) deallocate(mgc%pmat)
  if (allocated(mgc%cv))   deallocate(mgc%cv)

end subroutine free_mgchain


subroutine clhalf(n,mn,lmax,cl,bl,clh)
  implicit none
  integer, intent(in) :: n, mn, lmax
  double precision, intent(in) :: cl(n,0:lmax), bl(mn,0:lmax)
  double precision, intent(out) :: clh(n,mn,0:lmax,0:lmax)
  integer :: ni, l, mi

  clh = 0d0
  do ni = 1, n
    do l = 1, lmax
      do mi = 1, mn
        if (cl(ni,l)<=0d0) cycle
        clh(ni,mi,l,0:l) = dsqrt(cl(ni,l))*bl(mi,l)
      end do
    end do
  end do

end subroutine clhalf


subroutine correct_filtering(n,lmax,cl,filter,xlm)
  implicit none
  character(*), intent(in) :: filter
  integer, intent(in) :: n, lmax
  double precision, intent(in) :: cl(n,0:lmax)
  double complex, intent(inout) :: xlm(n,0:lmax,0:lmax)
  integer :: l, ni
  
  do l = 1, lmax
    do ni = 1, n
      if (cl(ni,l)<=0d0)  then
        xlm(ni,l,0:l) = 0d0
        cycle
      end if
      if (filter=='')   xlm(ni,l,0:l) = xlm(ni,l,0:l)/dsqrt(cl(ni,l))
      if (filter=='W')  xlm(ni,l,0:l) = xlm(ni,l,0:l)*dsqrt(cl(ni,l))
    end do
  end do

end subroutine


subroutine densemat_multi(n,lmax,arrn,r,imat,x)
  implicit none
  integer, intent(in) :: n, lmax, arrn
  double complex, intent(in) :: r(1:n,0:lmax,0:lmax)
  double complex, intent(in) :: imat(1:arrn,1:arrn)
  double complex, intent(out) :: x(1:n,0:lmax,0:lmax)
  integer :: i, j, l, m, p, q, ni, nj

  x = 0d0
  i = 0
  do ni = 1, n
    do l = 0, lmax
      do m = 0, l
        i = i + 1
        if (i>arrn) stop 'error'
        j = 0
        x(ni,l,m) = 0d0
        do nj = 1, n
          do p = 0, lmax
            do q = 0, p
              j = j + 1
              if (j>arrn) stop 'error'
              x(ni,l,m) = x(ni,l,m) + imat(i,j)*r(nj,p,q)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine densemat_multi


subroutine coarse_invmatrix(n,lmax,mgc)
  implicit none
  !I/O
  type(mg_chain), intent(inout) :: mgc
  integer, intent(in) :: n, lmax
  !double precision, intent(in), dimension(n,mn,0:lmax,0:lmax) :: clh
  !internal
  integer :: i, l, m, ni, mn
  double complex, allocatable :: x(:,:,:), Mx(:,:,:), A0(:,:), A1(:,:), y(:,:,:)

  mn = mgc%mnmax(mgc%n)
  mgc%matn = n*(lmax+1)*(lmax+2)/2

  ! prepare diagonal preconditioner
  allocate(mgc%imat(mgc%matn,mgc%matn),Mx(n,0:lmax,0:lmax),x(n,0:lmax,0:lmax))
  allocate(A0(mgc%matn,mgc%matn),A1(mgc%matn,mgc%matn)); A0=0d0; A1=0d0

  i = 0
  do ni = 1, n
    do l = 0, lmax
      do m = 0, l
        i = i + 1
        x = 0d0
        x(ni,l,m) = 1d0
        call matmul_lhs(n,mn,mgc%npix(mgc%n,:),lmax,mgc%cv(mgc%n,:),x,Mx)
        call trans_alm2array_1d(n,lmax,mgc%matn,Mx,A0(:,i))
        if (m/=0) then
          x = 0d0
          x(ni,l,m) = iu
          call matmul_lhs(n,mn,mgc%npix(mgc%n,:),lmax,mgc%cv(mgc%n,:),x,Mx)
          call trans_alm2array_1d(n,lmax,mgc%matn,Mx,A1(:,i))
          A0(:,i) = (A0(:,i)+A1(:,i)/iu)/2d0
        end if
      end do
    end do
  end do

  if (mgc%verbose) then
    i = 3*(3+1)/2
    write(*,*) '---'
    write(*,*) A0(i,i+1:i+5)
    write(*,*) '---'
    write(*,*) A0(i+1:i+5,i)
    write(*,*) '---'
    write(*,*) A0(i-1,i), A0(i,i-1)
    write(*,*) '---'
    write(*,*) A0(i-1,i+1), A0(i+1,i-1)
  end if

  call inv_matrix_c(A0,mgc%imat)

  !check
  if (mgc%verbose) then
    allocate(y(1:n,0:lmax,0:lmax)); y = 0d0
    y(1,3,2) = 1d0
    call densemat_multi(n,lmax,mgc%matn,y,mgc%imat,x)
    call matmul_lhs(n,mn,mgc%npix(mgc%n,:),lmax,mgc%cv(mgc%n,:),x,Mx)
    write(*,*) '---'
    write(*,*) Mx(1,3,0:4), Mx(1,4,0:4)
    deallocate(y)
    stop
  end if

  deallocate(Mx,x,A0,A1)


end subroutine coarse_invmatrix


subroutine cg_algorithm(n,lmax,b,x,mgc,chain,ratio)
!* Searching for a solution x of Ax = b with the Conjugate Gradient iteratively
!* The code assumes 
!*    1) A = [1 + C^1/2 N^-1 C^1/2]
!*    2) C^1/2 is diagonal
!*    3) N is diagonal in pixel space (statistically isotropic noise)
!*
  implicit none
  !I/O
  type(mg_chain), intent(inout) :: mgc
  integer, intent(in) :: n, lmax, chain
  !double precision, intent(in), dimension(n,mn,0:lmax,0:lmax) :: clh
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: b
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: x
  double precision, intent(out), optional :: ratio(mgc%itn(chain))
  !internal
  integer :: c, ni, mi, i, l
  double precision :: absb, absr, d, d0, td, alpha
  double precision :: mm(0:lmax,0:lmax)
  double complex, dimension(1:n,0:mgc%lmax(chain),0:mgc%lmax(chain)) :: r, z, p, Ap

  if (lmax/=mgc%lmax(chain)) stop 'error: lmax is inconsistent'
  absb = dsqrt(sum(abs(b)**2))


  if (chain==1) then
    !diagonal preconditioner
    write(*,*) 'precompute diagonal preconditioner'
    allocate(mgc%pmat(n,0:lmax,0:lmax))
    do ni = 1, n
      mm = 0d0
      do mi = 1, mgc%mnmax(1)
        !derived by averaging over m for int Ylm^* N^-1 Ylm
        mm = mm + mgc%cv(1,mi)%clh(ni,:,:)**2 * ave(mgc%cv(1,mi)%nij(ni,:)) 
      end do
      mgc%pmat(ni,:,:) = 1d0/(1d0+mm)
    end do
  end if

  !initial value (this is the solution if MA=I)
  call precondition(n,lmax,b,x,mgc,chain)

  !residual
  call matmul_lhs(n,mgc%mnmax(chain),mgc%npix(chain,:),lmax,mgc%cv(chain,:),x,r)
  r = b - r
  if (chain==1)  write(*,*) 'initial residual', dsqrt(sum(abs(r)**2)), '|b|^2', dsqrt(sum(abs(b)**2))
  !if (chain==1)  stop

  !set other values
  call precondition(n,lmax,r,z,mgc,chain)
  p = z

  !initial distance
  d0 = sum(conjg(r)*z)
  d  = d0

  do i = 1, mgc%itn(chain)

    call matmul_lhs(n,mgc%mnmax(chain),mgc%npix(chain,:),lmax,mgc%cv(chain,:),p,Ap)
    alpha = d/sum(conjg(p)*Ap)
    x = x + alpha*p
    r = r - alpha*Ap

    absr = dsqrt(sum(abs(r)**2))
    if (chain==1.and.mod(i,mgc%ro)==0)  write(*,*) i, absr/absb, d/d0
    if (chain==1.and.present(ratio))  ratio(i) = absr/absb

    call precondition(n,lmax,r,z,mgc,chain)
    td = sum(conjg(r)*z)
    p  = z + (td/d)*p
    d  = td

    ! check exit condition
    if (absr<mgc%eps(chain)*absb) then 
      !exit loop if |r|/|b| becomes very small
      if (chain==1) write(*,*) i, absr/absb
      exit
    end if

  end do

end subroutine cg_algorithm


subroutine precondition(n,lmax,r,x,mgc,chain)
  implicit none
  !I/O
  type(mg_chain), intent(inout) :: mgc
  integer, intent(in) :: n, lmax, chain
  !double precision, intent(in), dimension(n,mn,0:lmax,0:lmax) :: clh
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: r
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: x
  ! internal
  integer :: lmax0

  ! apply preconditioner
  if (mgc%n==1) then !diagonal preconditioning
    x = mgc%pmat(:,0:lmax,0:lmax)*r
  else
    !downgrade nij
    x = 0d0
    lmax0 = mgc%lmax(chain+1)
    if (chain+1==mgc%n) then ! multiply dense matrix at the coarsest gird
      call densemat_multi(n,lmax0,mgc%matn,r(:,0:lmax0,0:lmax0),mgc%imat,x(:,0:lmax0,0:lmax0))
    else
      call cg_algorithm(n,lmax0,r(:,0:lmax0,0:lmax0),x(:,0:lmax0,0:lmax0),mgc,chain+1) 
    end if
    if (lmax0+1<=mgc%lmax(chain))  x(:,lmax0+1:lmax,:) = mgc%pmat(:,lmax0+1:lmax,:)*r(:,lmax0+1:lmax,:)
  end if

end subroutine precondition


subroutine matmul_rhs(n,mn,npix,lmax,clh,cv,v)
  implicit none
  !I/O
  type(mg_covmat), intent(in) :: cv(mn)
  integer, intent(in) :: n, mn, npix(mn), lmax
  double precision, intent(in), dimension(n,mn,0:lmax,0:lmax) :: clh
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  integer :: mi
  double complex :: alm(n,0:lmax,0:lmax)

  v = 0d0
  do mi = 1, mn
    !nij = input map x N^-1
    call spht_map2alm(n,npix(mi),lmax,lmax,cv(mi)%nij,alm)
    v = v + clh(:,mi,:,:)*alm
  end do

end subroutine matmul_rhs


subroutine matmul_lhs(n,mn,npix2,lmax,cv,x,v)
  implicit none
  !I/O
  type(mg_covmat), intent(in) :: cv(mn)
  integer, intent(in) :: n, mn, npix2(mn), lmax
  !double precision, intent(in), dimension(n,mn,0:lmax,0:lmax) :: clh
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: x
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  integer :: mi, npix
  double complex :: alm(n,0:lmax,0:lmax)

  !make 1+C^{1/2}N^{-1}C^{1/2}
  v = x
  do mi = 1, mn
    npix = size(cv(mi)%nij,dim=2)
    call matmul_cninv(n,npix,lmax,cv(mi)%clh,cv(mi)%nij,x,alm)
    v = v + alm
  end do

end subroutine matmul_lhs


subroutine matmul_cninv(n,npix,lmax,clh,nij,x,v)
  implicit none
  !I/O
  integer, intent(in) :: n, npix, lmax
  double precision, intent(in), dimension(n,0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: x
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  double precision :: map(n,0:npix-1)
  double complex :: alm(n,0:lmax,0:lmax)

  alm = clh*x
  !multiply noise covariance in pixel space
  call spht_alm2map(n,npix,lmax,lmax,alm,map)
  map = map*nij
  call spht_map2alm(n,npix,lmax,lmax,map,alm)
  !multiply C^{1/2}
  v = clh*alm

end subroutine matmul_cninv


end module hp_cgd


