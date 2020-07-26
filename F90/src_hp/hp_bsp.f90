!////////////////////////////////////////////////////!
! Bispectrum calcuation in fullsky
!////////////////////////////////////////////////////!

module hp_bsp
  use alm_tools, only: alm2map, map2alm
  use constants, only: pi

  private alm2map
  private pi

contains 


subroutine equi(lmin,lmax,alm,bispec,bst)
! Compute equilateral shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l,l]
  implicit none
  !I/O
  integer, intent(in) :: lmin, lmax, bst
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:)
  double complex, allocatable :: klm(:,:,:)

  nside = bst*2**(int(dlog(dble(lmax))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1),klm(1,0:lmax,0:lmax))

  klm = 0d0
  do l = lmin, lmax !multipole filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,lmax,lmax,klm,kmap)
  bispec = sum(kmap**3)*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine equi


subroutine fold(lmin,lmax,alm,bispec,bst)
! Compute folded shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l/2,l/2]
  implicit none
  !I/O
  integer, intent(in) :: lmin, lmax, bst
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nside = bst*2**(int(dlog(dble(lmax))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:lmax,0:lmax))

  klm = 0d0
  do l = lmin, lmax !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = max(2,int(lmin/2d0)), int(lmax/2d0)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,lmax,lmax,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,lmax,lmax,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine fold


subroutine sque(eL,sL,l1,alm,bispec,bst)
! Compute squeezed shape of the unnormalized binned reduced bispectrum for a given alm, b[sL,eL,eL]
  implicit none
  !I/O
  integer, intent(in) :: l1, bst
  integer, intent(in), dimension(2) :: eL, sL
  double complex, intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  if (max(sL(2),eL(2))>l1) stop 'error (sque): l1 is too small'

  nside = bst*2**(int(dlog(dble(l1))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:l1,0:l1))

  klm = 0d0
  do l = sL(1), sL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = eL(1), eL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,l1,l1,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,l1,l1,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine sque


subroutine isos(eL,aL,l1,alm,bispec,bst)
! Compute isosceles shape of the unnormalized binned reduced bispectrum for a given alm, b[eL,aL,aL]
  implicit none
  !I/O
  integer, intent(in) :: l1, bst
  integer, intent(in), dimension(2) :: eL, aL
  double complex, intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  if (max(aL(2),eL(2))>l1) stop 'error (isos): l1 is too small'

  nside = bst*2**(int(dlog(dble(l1))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:l1,0:l1))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = aL(1), aL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,l1,l1,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,l1,l1,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine isos


subroutine xequi(lmin,lmax,n,alm,bispec,bst)
! Compute equilateral shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l,l]
  implicit none
  !I/O
  integer, intent(in) :: lmin, lmax, bst, n
  double complex, intent(in), dimension(1:n,0:lmax,0:lmax) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, i, nside
  double precision :: bsp
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nside = bst*2**(int(dlog(dble(lmax))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,n),klm(n,0:lmax,0:lmax))

  klm = 0d0
  do l = lmin, lmax !multipole filtering
    klm(:,l,0:l) = alm(:,l,0:l)
  end do
  do i = 1, n
    call alm2map(nside,lmax,lmax,klm(i:i,:,:),kmap(:,i))
  end do

  select case(n)
  case(1)
    bsp = sum(kmap(:,1)**3)
  case(2)
    bsp = sum(kmap(:,1)*kmap(:,2)**2)
  case(3)
    bsp = sum(kmap(:,1)*kmap(:,2)*kmap(:,3))
  case default
    stop 'error (xequi): n should be 1<=n<=3'
  end select

  bispec = bsp*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine xequi


subroutine xfold(lmin,lmax,n,alm,bispec,bst)
! Compute folded shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l/2,l/2]
  implicit none
  !I/O
  integer, intent(in) :: lmin, lmax, bst, n
  double complex, intent(in), dimension(1:n,0:lmax,0:lmax) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, i, nn, nside
  double precision :: bsp
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nn = max(n,2)
  nside = bst*2**(int(dlog(dble(lmax))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,nn),klm(nn,0:lmax,0:lmax))

  klm = 0d0
  do l = lmin, lmax !ell filtering
    klm(1,l,0:l) = alm(1,l,0:l)
  end do
  do l = max(2,int(lmin/2d0)), int(lmax/2d0)
    if (n==1)  klm(2,l,0:l)   = alm(1,l,0:l)
    if (n>=2)  klm(2:n,l,0:l) = alm(2:n,l,0:l)
  end do

  do i = 1, nn
    call alm2map(nside,lmax,lmax,klm(i:i,:,:),kmap(:,i))
  end do

  select case(n)
  case(1)
    bsp = sum(kmap(:,1)**3)
  case(2)
    bsp = sum(kmap(:,1)*kmap(:,2)**2)
  case(3)
    bsp = sum(kmap(:,1)*kmap(:,2)*kmap(:,3))
  case default
    stop 'error (xfold): n should be 1<=n<=3'
  end select

  bispec = bsp*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine xfold


subroutine xsque(eL,sL,l1,n,alm,bispec,bst)
! Compute squeezed shape of the unnormalized binned reduced bispectrum for a given alm, b[sL,eL,eL]
  implicit none
  !I/O
  integer, intent(in) :: l1, bst, n
  integer, intent(in), dimension(2) :: eL, sL
  double complex, intent(in), dimension(1:n,0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, i, nside, nn
  double precision :: bsp
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nn = max(n,2)
  nside = bst*2**(int(dlog(dble(l1))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,nn),klm(nn,0:l1,0:l1))

  klm = 0d0
  do l = sL(1), sL(2) !ell filtering
    klm(1,l,0:l) = alm(1,l,0:l)
  end do
  do l = eL(1), eL(2)
    if (n==1)  klm(2,l,0:l)   = alm(1,l,0:l)
    if (n>=2)  klm(2:n,l,0:l) = alm(2:n,l,0:l)
  end do

  do i = 1, nn
    call alm2map(nside,l1,l1,klm(i:i,:,:),kmap(:,i))
  end do

  select case(n)
  case(1)
    bsp = sum(kmap(:,1)**3)
  case(2)
    bsp = sum(kmap(:,1)*kmap(:,2)**2)
  case(3)
    bsp = sum(kmap(:,1)*kmap(:,2)*kmap(:,3))
  case default
    stop 'error (xsque): n should be 1<=n<=3'
  end select

  bispec = bsp*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine xsque


subroutine xisos(eL,aL,l1,n,alm,bispec,bst)
! Compute isosceles shape of the unnormalized binned reduced bispectrum for a given alm, b[eL,aL,aL]
  implicit none
  !I/O
  integer, intent(in) :: l1, bst, n
  integer, intent(in), dimension(2) :: eL, aL
  double complex, intent(in), dimension(1:n,0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, i, nside, nn
  double precision :: bsp
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  if (max(aL(2),eL(2))>l1) stop 'error (isos): l1 is too small'

  nn = max(n,2)
  nside = bst*2**(int(dlog(dble(l1))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1,nn),klm(nn,0:l1,0:l1))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(1,l,0:l)
  end do
  do l = aL(1), aL(2)
    if (n==1)  klm(2,l,0:l)   = alm(1,l,0:l)
    if (n>=2)  klm(2:n,l,0:l) = alm(2:n,l,0:l)
  end do

  do i = 1, nn
    call alm2map(nside,l1,l1,klm(i:i,:,:),kmap(:,i))
  end do

  select case(n)
  case(1)
    bsp = sum(kmap(:,1)**3)
  case(2)
    bsp = sum(kmap(:,1)*kmap(:,2)**2)
  case(3)
    bsp = sum(kmap(:,1)*kmap(:,2)*kmap(:,3))
  case default
    stop 'error (xisos): n should be 1<=n<=3'
  end select

  bispec = bsp*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine xisos


end module hp_bsp


