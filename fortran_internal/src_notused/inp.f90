!//////////////////////////////////////////////////////////////////!
! * Conjugate Gradient Algorithm for Gaussian-Constrained Inpainting
! * Not well tested
!//////////////////////////////////////////////////////////////////!

module mycinv
  use myconst, only: ac2rad
  use myfftw, only: dft_pol
  implicit none

  !* local parameter
  integer, parameter :: dlc = KIND(0d0)

  private dft_pol, dlc, ac2rad

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
  complex(dlc), intent(in) :: alm(:,:)
! [outputs]
!   ilm --- inverse variance filteres multipole
  complex(dlc), intent(out) :: ilm(:,:)
! [internal]
  integer :: npix
  complex(dlc), allocatable :: xlm(:,:), ylm(:,:)

  npix = nn(1)*nn(2)
  allocate(xlm(2,npix),ylm(2,npix)); xlm=0d0; ylm=0.d0
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
  complex(dlc), intent(in) :: alm(:,:)
  complex(dlc), intent(out) :: xlm(:,:)
  !internal
  integer :: i, npix
  double precision, allocatable :: CF(:,:)
  complex(dlc), allocatable :: b(:,:)

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
  complex(dlc), intent(in) :: b(:,:)
  complex(dlc), intent(out) :: x(:,:)
  !internal
  integer :: i, j, iter, roundoff=50, simn, npix
  double precision :: dr, dr0, tdr, alpha, eps = 1d-6, sigp(2)
  complex(dlc), allocatable, dimension(:,:) :: r, tAr, p, r0, Ap

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
  complex(dlc), intent(in) :: x(:,:)
  complex(dlc), intent(out) :: v(:,:)
  !internal
  integer :: npix
  double precision, allocatable :: QU(:,:)
  complex(dlc), allocatable :: EB(:,:)

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
! * compute right-hand-side of Ax=b
!     b = C^1/2 N^-1 x
  implicit none
  !I/O
  integer, intent(in) :: nn(:)
  double precision, intent(in) :: D(2), C(:,:), N(:,:)
  complex(dlc), intent(in) :: dat(:,:)
  complex(dlc), intent(out) :: b(:,:)
  !internal
  integer :: npix
  double precision, allocatable :: QU(:,:)
  complex(dlc), allocatable :: EB(:,:)

  npix = nn(1)*nn(2)
  allocate(EB(2,npix),QU(npix,2))
  EB = dat
  call dft_pol(QU,nn,D,EB,-1)
  QU = QU*N
  call dft_pol(QU,nn,D,EB,1)
  b = C*EB
  deallocate(EB,QU)

end subroutine vecb

end module mycinv


#ifdef useold
!//////////////////////////////////////////////////////////////////!
! * Old subroutines
!//////////////////////////////////////////////////////////////////!

!this routine work only for one frequency case (nc=1)
subroutine inpaint_interface(nn,D,CL2D,NTP,mask,T,P,f1,f2,DoInp)
  implicit none
  !I/O
  logical, intent(in), optional :: DoInp
  character(*), intent(in), optional :: f1, f2
  integer, intent(in) :: nn(2)
  real(dl), intent(in) :: D(2), CL2D(:,:), NTP(:,:)
  complex(dlc), intent(in), optional :: T(:,:), P(:,:)
  complex(dlc), intent(in) :: mask(:)
  !internal
  integer :: i, n, npix
  complex(dlc), allocatable :: xT(:), xP(:), yT(:), yP(:), inpT(:), inpP(:)

  npix = size(Cl2D,dim=2)

  allocate(xT(npix),xP(npix),yT(npix),yP(npix),inpT(npix),inpP(npix)); yT=0d0; yP=0d0
  if(present(T)) then
    if(present(P)) then
      call EstField(nn,D,CL2D,NTP,mask,xT,xP,T,P)
      if(present(DoInp)) call EstField(nn,D,CL2D,NTP,mask,yT,yP,T,P,.true.)
    else
      call EstField(nn,D,CL2D,NTP,mask,xT,xP,T)
      if(present(DoInp)) call EstField(nn,D,CL2D,NTP,mask,yT,yP,T,fluc=.true.)
    end if
  else if(present(P)) then
    call EstField(nn,D,CL2D,NTP,mask,xT,xP,P)
    if(present(DoInp)) call EstField(nn,D,CL2D,NTP,mask,yT,yP,P,fluc=.true.)
  end if

  inpT = xT+yT
  inpP = xP+yP
  !if(present(f1).and.present(f2)) call InpaintedMap(nn,D,inpT,1,f1,f2)

  deallocate(xT,yT,xP,yP,inpT,inpP)

end subroutine inpaint_interface


subroutine EstField(nn,D,CL2D,NTP,mask,xT,xP,T,P,fluc)
  implicit none
  !I/O
  logical, intent(in), optional :: fluc
  integer, intent(in) :: nn(2)
  real(dl), intent(in) :: D(2), CL2D(:,:), NTP(:,:)
  complex(dlc), intent(in), optional :: T(:,:), P(:,:)
  complex(dlc), intent(out) :: xT(:), xP(:)
  complex(dlc), intent(in) :: mask(:)
  !internal
  integer :: i, n, nu, nc, l, num, npix
  complex(dlc) :: E, B
  complex(dlc), allocatable :: bT(:), bP(:), T0(:), P0(:), T1(:,:), P1(:,:)

  nc   = size(NTP)
  npix = size(Cl2D,dim=2)
  allocate(bT(npix),bP(npix),T0(npix),P0(npix),T1(nc,npix),P1(nc,npix))

  if(present(fluc)) then !fluctuations
    call InitRandom(-1)
    do n = 1, npix
      if(present(T)) T0(n) = Gaussian1()
      if(present(P)) P0(n) = Gaussian1()+iu*Gaussian1() !need check
      do nu = 1, nc
        if(present(T)) T1(nu,n) = Gaussian1()
        if(present(P)) P1(nu,n) = Gaussian1()+iu*Gaussian1() !need check
      end do
    end do
  else !mean field
    if(present(T)) T0 = 0d0
    if(present(P)) P0 = 0d0
    if(present(T)) T1 = T
    if(present(P)) P1 = P
  end if

  do n = 1, npix
    do nu = 1, nc
      if(present(T)) xT(n) = T1(nu,n)*dsqrt(Cl2D(1,n))/(Cl2D(1,n)+NTP(1,nu)**2)
      if(present(P)) then
        E = (P1(nu,n)+conjg(P1(nu,n)))*0.5d0
        B = (P1(nu,n)-conjg(P1(nu,n)))*0.5d0/iu
        xP(n) = E*dsqrt(Cl2D(2,n))/(Cl2D(2,n)+NTP(2,nu)**2) + iu*B*dsqrt(Cl2D(3,n))/(Cl2D(3,n)+NTP(2,nu)**2)
      end if
    end do
  end do

  call vecb(nn,D,bT,T0,T1,fluc,mask,CL2D,NTP)
  call cg_algorithm(nn,D,bT,xT,mask,num,CL2D,NTP)

  if(if(present(T))) xT = xT*dsqrt(Cl2D(1,:))
  if(if(present(P))) xP = (xP+conjg(xP))*dsqrt(Cl2D(2,:))*0.5d0 + iu*(xP-conjg(xP))*dsqrt(Cl2D(3,:))*0.5d0/iu

end subroutine EstField 


subroutine cg_algorithm(nn,D,b,x,beam,iNl,nstp,eps)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  integer, intent(in), optional :: nstp
  real(dl), intent(in) :: D(2), beam(:,:,:), iNl(:,:,:)
  real(dl), intent(in), optional :: eps
  complex(dlc), intent(in) :: b(:,:)
  complex(dlc), intent(out) :: x(:,:)
  !internal
  character(LEN=50) :: iterc, f
  logical :: diag
  integer :: i, j, iter, roundoff=50, simn, nc, npix
  real(dl) :: eps0, delta, tdelta, alpha
  complex(dlc), allocatable, dimension(:,:) :: r, tAr, p, r0, Ap, map

  num  = size(x,dim=1)
  npix = nn(1)*nn(2)
  allocate(r(num,npix),tAr(num,npix),p(num,npix),r0(num,npix),Ap(num,npix),map(num,npix))

  nc   = 1
  diag = .true.

  call asub(nn,D,x,r,beam,iNl)
  r = b - r
  call pre_op(r,p,beam,iNl)
  delta = sum(conjg(r)*p)

  iter = 0
  simn = 1000
  eps0 = 1d-5
  if(present(nstp)) simn = nstp
  if(present(eps))  eps0 = eps

  do i = 1,simn

    call asub(nn,D,p,Ap,beam,iNl)
    alpha = delta/sum(conjg(p)*Ap)

    x = x + alpha*p

    iter = iter + 1
    if(iter-int(dble(iter)/dble(roundoff))*roundoff==0) then 
      r0 = r
      call asub(nn,D,x,r0,beam,iNl)
      r0 = b - r0
      !write(*,*) sum(conjg(r0(1,:))*r0(1,:))/dble(nn(1)*nn(2))
      !write(*,*) sum(conjg(r0(2,:))*r0(2,:))/dble(nn(1)*nn(2))
    end if
    r = r - alpha*Ap

    call pre_op(r,tAr,beam,iNl)
    tdelta = sum(conjg(r)*tAr)

    p = tAr + (tdelta/delta)*p
    delta = tdelta
    
    if(delta<eps0**2) exit !Norm of r is sufficiently small
    write(*,*) iter, delta
 
  end do

  deallocate(r, tAr, p, r0, Ap, map)

end subroutine cg_algorithm


subroutine asub(nn,D,x,v,beam,iNl)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  real(dl), intent(in) :: D(2), beam(:,:,:), iNl(:,:,:)
  complex(dlc), intent(in) :: x(:,:)
  complex(dlc), intent(out) :: v(:,:)
  !internal
  integer :: nu
  complex(dlc), allocatable :: cmb(:,:)

  allocate(cmb(3,nn(1)*nn(2)))
  v = 0d0
  do nu = 1, size(beam,dim=1)
    cmb = beam(nu,:,:)*x
    call dft_all(cmb,nn,D,-1)
    cmb = cmb*iNl(nu,:,:)
    call dft_all(cmb,nn,D,1)
    v = v + beam(nu,:,:)*cmb
  end do
  v = v + x
  deallocate(cmb)

end subroutine asub


subroutine vecb(nn,D,b,v0,v1,beam,iNl)
  implicit none
  !I/O
  integer, intent(in) :: nn(:)
  real(dl), intent(in) :: D(2), beam(:,:,:), iNl(:,:,:)
  complex(dl), intent(out) :: b(:,:)
  complex(dlc), intent(in) :: v0(:,:), v1(:,:,:)
  !internal
  integer :: nu
  complex(dlc), allocatable :: cmb(:)

  allocate(cmb(3,nn(1)*nn(2)))
  b = 0d0
  do nu = 1, size(v1,dim=2)
    cmb = v1(nu,:,:)
    call dft_all(cmb,nn,D,-1)
    cmb = cmb*iNl(nu,:,:)
    call dft_all(cmb,nn,D,1)
    b = b + beam(nu,:,:)*cmb
  end do
  b = b + v0
  deallocate(cmb)

end subroutine vecb


subroutine pre_op(nn,D,x,v,beam,iNl)
  implicit none 
  !I/O
  real(dl), intent(in) :: beam(:,:,:), iNl(:,:,:)
  complex(dlc), intent(in) :: x(:,:)
  complex(dlc), intent(out) :: v(:,:)
  !internal
  integer :: nu, nc

  nc = 1
  !Multiply preconditioner
  do nu = 1, nc
    v = x/(iNl(nu,:,:)*beam(nu,:,:) + 1d0)
  end do

end subroutine pre_op


subroutine InpaintedMap(nn,D,inp,res,f1,f2)
  implicit none 
  !I/O
  character(*), intent(in) :: f1, f2
  integer, intent(in) :: res, nn(2)
  real(dl), intent(in) :: D(2)
  complex(dlc), intent(in) :: inp(:,:)
  !internal
  integer :: i, l
  complex(dlc) :: Xp(nn(1)*nn(2)),Xm(nn(1)*nn(2))
  complex(dlc) :: MAP(2,nn(1),nn(2)),TQU(3,nn(1),nn(2))

  Xp(:) = inp(1,:) + iu*inp(2,:)
  Xm(:) = inp(1,:) - iu*inp(2,:)
  ! Fourier -> Real 
  call dft_pol(Xp,nn,D,-1)
  call dft_pol(Xm,nn,D,-1)
  call ArrayToMap(nn,Xp,MAP(1,:,:))
  call ArrayToMap(nn,Xm,MAP(2,:,:))

  TQU(1,:,:) = 0d0
  TQU(2,:,:) = (MAP(1,:,:)+MAP(2,:,:))*cmplx(Tcmb)*0.5d0
  TQU(3,:,:) = (MAP(1,:,:)-MAP(2,:,:))*cmplx(Tcmb)*0.5d0/iu
  call Map_Write(nn,res,f1,MAP2D=dble(TQU(1,:,:)))

end subroutine InpaintedMap
#endif useold

