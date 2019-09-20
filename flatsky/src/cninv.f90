!//////////////////////////////////////////////////////////////////!
! Conjugate Gradient Decent Algorithm
!//////////////////////////////////////////////////////////////////!

module cninv
  !use myconst, only: ac2rad
  !use myfftw, only: dft_pol
  implicit none

  !private dft_pol, dlc, ac2rad

contains 

subroutine cinv_pol(nn,D,alm,C,N,F,ilm)
! * computing inverse-variance filtered multipoles: C^-1d
  implicit none
! [inputs]  
!   nn (2)      --- number of x and y grids
!   D(2)        --- each side length
!   alm(2,npix) --- observed Fourier mode, assuming the form, alm = F*cmb + noise
!   C(2,npix)   --- half of the signal covariance in Fourier space: C^1/2
!   N(npix,2)   --- noise covariance in real space
!   F(2,npix)   --- filter operation in Fourier space
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), C(:,:), N(:,:), F(:,:)
  double complex, intent(in) :: alm(:,:)
! [outputs]
!   ilm --- inverse variance filteres multipole
  double complex, intent(out) :: ilm(:,:)
! [internal]
  integer :: npix
  double complex, allocatable :: xlm(:,:)

  npix = nn(1)*nn(2)
  allocate(xlm(2,npix)); xlm=0d0
  if(size(N,dim=2)>2) stop "error: dimension of invN is stragne"
  if(size(C,dim=1)>2) stop "error: dimension of cov is stragne"
  call est_mean(nn,D,alm,C,N,F,xlm)
  ilm = xlm+ylm
  deallocate(xlm,ylm)

end subroutine cinv_pol


subroutine est_mean(nn,D,alm,C,N,F,xlm)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), C(:,:), N(:,:), F(:,:)
  double complex, intent(in) :: alm(:,:)
  double complex, intent(out) :: xlm(:,:)
  !internal
  integer :: i, npix
  double precision, allocatable :: CF(:,:)
  double complex, allocatable :: b(:,:)

  npix = nn(1)*nn(2)
  allocate(b(2,npix),CF(2,npix));  b=0d0;  CF=0d0
  CF = C*F ! CF = C^1/2 * F
  call vecb(nn,D,alm,CF,N,b)
  call CG_algorithm(nn,D,b,CF,N,xlm)
  do i = 1, npix
    if(C(1,i)>0d0) xlm(:,i) = xlm(:,i)/C(:,i)
  end do
  deallocate(b,CF)

end subroutine est_mean


subroutine CG_algorithm(nn,D,b,C,N,x,forcestop)
! * Searching a solution of Ax=b with the Conjugate Gradient Algorithm iteratively, where A is almost diagonal
  implicit none
  !I/O
  integer, intent(in) :: nn(:)
  integer, intent(in), optional :: forcestop
  double precision, intent(in) :: D(2), C(:,:), N(:,:)
  double complex, intent(in) :: b(:,:)
  double complex, intent(out) :: x(:,:)
  !internal
  integer :: i, j, iter, roundoff=50, simn, npix
  double precision :: dr, dr0, tdr, alpha, eps = 1d-6, sigp(2)
  double complex, allocatable, dimension(:,:) :: r, tAr, p, r0, Ap

  npix = nn(1)*nn(2)
  sigp = sum(N,dim=1)/dble(npix)
  allocate(r(2,npix),tAr(2,npix),p(2,npix),r0(2,npix),Ap(2,npix))

  call mat_multiplication(nn,D,x,C,N,r)  ! r = A*x
  r = b - r
  !multiply preconditioner
  do j = 1, npix
    p(:,j) = r(:,j)/(1d0+sigp(:)*C(:,j)**2)   ! p = tA*r
  end do

  dr0 = sum(conjg(r)*p)
  dr = dr0

  iter = 0
  simn = 10000
  if(present(forcestop)) simn = forcestop

  do i = 1, simn

    call mat_multiplication(nn,D,p,C,N,Ap)  ! Ap = A*p
    alpha = dr/sum(conjg(p)*Ap)
    x = x + alpha*p
    iter = iter + 1

    if(iter-int(dble(iter)/dble(roundoff))*roundoff==0) then 
      r0 = r
      call mat_multiplication(nn,D,x,C,N,r0)
      r0 = b - r0
      write(*,*) sum(conjg(r0(1,:))*r0(1,:))/dble(npix), sum(conjg(r0(2,:))*r0(2,:))/dble(npix)
    end if

    r = r - alpha*Ap
    !* tAr = tA*r
    do j = 1, npix
      tAr(:,j) = r(:,j)/(1d0+sigp(:)*C(:,j)**2)
    end do
    tdr = sum(conjg(r)*tAr)
    p = tAr + (tdr/dr)*p
    dr = tdr

    if(dr<eps**2*dr0) then 
      write(*,*) iter, dr/dr0
      exit !Norm of r is sufficiently small
    end if 

  end do

  deallocate(r, tAr, p, r0, Ap)

end subroutine CG_algorithm


subroutine mat_multiplication(nn,D,x,C,N,v)
! * multiplying matrix
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), C(:,:), N(:,:)
  double complex, intent(in) :: x(:,:)
  double complex, intent(out) :: v(:,:)
  !internal
  integer :: npix
  double precision, allocatable :: QU(:,:)
  double complex, allocatable :: EB(:,:)

  npix = nn(1)*nn(2)
  allocate(EB(2,npix),QU(npix,2))
  EB = C*x
  call dft_pol(QU,nn,D,EB,-1)
  QU = QU*N
  call dft_pol(QU,nn,D,EB,1)
  v = x + C*EB
  deallocate(EB,QU)

end subroutine mat_multiplication


subroutine vecb(nn,D,dat,C,N,b)
! compute R.H.S. of Ax = b where b = C^1/2 N^-1 x
  implicit none
  !I/O
  integer, intent(in) :: nn(:)
  double precision, intent(in) :: D(2), C(:,:), N(:,:)
  double complex, intent(in) :: dat(:,:)
  double complex, intent(out) :: b(:,:)
  !internal
  integer :: npix
  double precision, allocatable :: QU(:,:)
  double complex, allocatable :: EB(:,:)

  npix = nn(1)*nn(2)
  allocate(EB(2,npix),QU(npix,2))
  EB = dat
  call dft_pol(QU,nn,D,EB,-1)
  QU = QU*N
  call dft_pol(QU,nn,D,EB,1)
  b = C*EB
  deallocate(EB,QU)

end subroutine vecb

end module cninv


