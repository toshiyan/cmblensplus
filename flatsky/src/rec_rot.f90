!///////////////////////////////////////////////////////////!
! * Rotation Reconstruction Kernel
!///////////////////////////////////////////////////////////!

module rotrecflat
  use myconst, only: iu, i1, pi, dlc, twopi
  use anaflat, only: elarrays
  use myfftw,  only: dft
  implicit none

  private iu, i1, pi, dlc, twopi
  private elarrays
  private dft

contains 


!/////////////////////////////////////////////////////////////////////////////////////!
! cosmic birefringence estimator
!

subroutine quadeb_cb(nn,D,E,B,EE,BB,rL,olm,eL)
  implicit none
! [input]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   EE(:)   --- E-mode Cl in 2D grids
!   BB(:)   --- B-mode Cl in 2D grids
!   E(:)    --- E-modes
!   B(:)    --- B-modes
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), BB(:), D(1:2)
  complex(dlc), intent(in) :: E(:), B(:)
!
! (optional)
!   eL(1:2) --- multipole range of output
!   olm(:)  --- rotation
  integer, intent(in), optional :: eL(1:2)
  complex(dlc), intent(out), optional :: olm(:)
!
! [internal]
  integer :: n, npix
  double precision :: els(nn(1)*nn(2))
  complex(dlc) :: ei2p(nn(1)*nn(2))
  complex(dlc), allocatable :: wE(:), wB(:), alm(:), blm(:)

  npix = nn(1)*nn(2)
  call elarrays(nn,D,els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (quadeb_cb): size of EE /= npix.'

  !* filtering
  allocate(wE(npix),wB(npix)); wE=0d0; wB=0d0
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wE(n) = 2d0*EE(n)*E(n)*ei2p(n)
    wB(n) = B(n)*conjg(ei2p(n))
  end do 

  !* convolution
  call dft(wE,nn,D,-1)
  call dft(wB,nn,D,-1)
  allocate(alm(npix))
  alm = dble(wE*wB)
  deallocate(wE,wB)
  call dft(alm,nn,D,1)

  allocate(blm(npix),wE(npix),wB(npix)); wE=0d0; wB=0d0; blm=0d0
  if (size(BB)/=npix)  stop 'error (quadeb_cb): size of BB /= npix.'
  if (sum(abs(BB))==0d0) goto 1
  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wE(n) = E(n)*ei2p(n)
    wB(n) = 2d0*BB(n)*B(n)*conjg(ei2p(n))
  end do 
  !* convolution
  call dft(wE,nn,D,-1)
  call dft(wB,nn,D,-1)
  blm = dble(wE*wB)
  call dft(blm,nn,D,1)
1 deallocate(wE,wB)

  !* estimator
  do n = 1, npix
    if(present(eL).and.(els(n)<eL(1).or.els(n)>eL(2))) alm(n)=0d0
    if(present(eL).and.(els(n)<eL(1).or.els(n)>eL(2))) blm(n)=0d0
    if(present(olm)) olm(n) = alm(n) - blm(n)
  end do
  deallocate(alm,blm)

end subroutine quadeb_cb


subroutine quadeb_cb_nobb(nn,D,E,B,EE,rL,olm,eL)
  implicit none
! [input]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   fC(:)   --- E-mode Cl in 2D grids
!   E(:)    --- E-modes
!   B(:)    --- B-modes
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), D(1:2)
  complex(dlc), intent(in) :: E(:), B(:)
!
! (optional)
!   eL(1:2) --- multipole range of output
!   olm(:)  --- rotation
  integer, intent(in), optional :: eL(1:2)
  complex(dlc), intent(out), optional :: olm(:)
!
! [internal]
  integer :: n, npix
  double precision :: els(nn(1)*nn(2))
  complex(dlc) :: ei2p(nn(1)*nn(2))
  complex(dlc), allocatable :: wE(:), wB(:), alm(:)

  npix = nn(1)*nn(2)
  call elarrays(nn,D,els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (quadeb_cb): size of EE /= npix.'

  !* filtering
  allocate(wE(npix),wB(npix)); wE=0d0; wB=0d0
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wE(n) = 2d0*EE(n)*E(n)*ei2p(n)
    wB(n) = B(n)*conjg(ei2p(n))
  end do 

  !* convolution
  call dft(wE,nn,D,-1)
  call dft(wB,nn,D,-1)
  allocate(alm(npix))
  alm = dble(wE*wB)
  deallocate(wE,wB)
  call dft(alm,nn,D,1)

  !* estimator 
  do n = 1, npix
    if(present(eL).and.(els(n)<eL(1).or.els(n)>eL(2))) alm(n)=0d0
    if(present(olm)) olm(n) = alm(n)
  end do
  deallocate(alm)

end subroutine quadeb_cb_nobb


subroutine alflat_eb_cb(nn,D,IE,IB,EE,BB,rL,Ala,eL)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), rL(1:2)
  integer, intent(in), optional :: eL(1:2)
  double precision, intent(in) :: EE(:), BB(:), IE(:), IB(:), D(1:2)
  double precision, intent(out), optional :: Ala(:)
! [internal]
  integer :: i, n, npix
  double precision :: els(nn(1)*nn(2)), iAl
  complex(dlc), allocatable :: Al(:,:), Bl(:,:), Cl(:,:), A1(:,:), A2(:,:), ei2p(:)

  write(*,*) 'EB flat (CB)'

  npix = nn(1)*nn(2)
  allocate(ei2p(npix))
  call elarrays(nn,D,els=els,ei2p=ei2p)
  if (size(EE)/=npix)  stop 'error (alflat_eb_cb): size(EE) /= npix.'
  if (size(BB)/=npix)  stop 'error (alflat_eb_cb): size(BB) /= npix.'

  !* EE^2 part
  allocate(A1(2,npix),A2(2,npix));  A1=0d0;  A2=0d0
  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(:,n) = IE(n)*EE(n)**2 * [i1,ei2p(n)**2]
    A2(:,n) = IB(n) * [i1,conjg(ei2p(n))**2]
  end do
  !* convolution
  allocate(Al(2,npix))
  do i = 1, 2
    call dft(A1(i,:),nn,D,-1)
    call dft(A2(i,:),nn,D,-1)
    Al(i,:) = A2(i,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
    Al(i,:) = dble(Al(i,:))
  end do
  deallocate(A1,A2)

  !* BB^2 part
  allocate(Bl(2,npix),A1(2,npix),A2(2,npix));  A1=0d0;  A2=0d0;  Bl=0d0
  if (sum(abs(BB))==0d0) goto 1
  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(:,n) = IE(n) * [i1,ei2p(n)**2]
    A2(:,n) = IB(n)*BB(n)**2 * [i1,conjg(ei2p(n))**2]
  end do
  !* convolution
  do i = 1, 2
    call dft(A1(i,:),nn,D,-1)
    call dft(A2(i,:),nn,D,-1)
    Bl(i,:) = A2(i,:)*A1(i,:)
    call dft(Bl(i,:),nn,D,1)
    Bl(i,:) = dble(Bl(i,:))
  end do
1 deallocate(A1,A2)


  !* EExBB part
  allocate(Cl(2,npix),A1(2,npix),A2(2,npix));  A1=0d0;  A2=0d0;  Cl=0d0
  if (sum(abs(BB))==0d0) goto 2
  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(:,n) = EE(n)*IE(n) * [i1,ei2p(n)**2]
    A2(:,n) = BB(n)*IB(n) * [i1,conjg(ei2p(n))**2]
  end do
  !* convolution
  do i = 1, 2
    call dft(A1(i,:),nn,D,-1)
    call dft(A2(i,:),nn,D,-1)
    Cl(i,:) = A2(i,:)*A1(i,:)
    call dft(Cl(i,:),nn,D,1)
    Cl(i,:) = dble(Cl(i,:))
  end do
2 deallocate(A1,A2)

  !* normalization
  do n = 1, npix
    if (present(eL).and.els(n)<eL(1).or.els(n)>eL(2))  cycle
    iAl = 2d0*( sum(Al(1:2,n)) + sum(Bl(1:2,n)) ) - 4d0*sum(Cl(1:2,n))
    if (iAl<=0d0) cycle
    if (present(Ala))  Ala(n) = 1d0 / iAl
  end do

  deallocate(Al,Bl,Cl,ei2p)

end subroutine alflat_eb_cb


subroutine alflat_eb_cb_nobb(nn,D,IE,IB,EE,rL,Ala,eL)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), rL(1:2)
  integer, intent(in), optional :: eL(1:2)
  double precision, intent(in) :: EE(:), IE(:), IB(:), D(1:2)
  double precision, intent(out), optional :: Ala(:)
! [internal]
  integer :: i, n, npix
  double precision :: els(nn(1)*nn(2)), iAl
  complex(dlc), allocatable :: Al(:,:), A1(:,:), A2(:,:), ei2p(:)

  write(*,*) 'EB flat (CB)'

  npix = nn(1)*nn(2)
  allocate(A1(2,npix),A2(2,npix),ei2p(npix));  A1=0d0;  A2=0d0
  call elarrays(nn,D,els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (alflat_eb_cb): size of EE is strange.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(:,n) = [ei2p(n)**2,i1] * EE(n)**2 * IE(n)
    A2(:,n) = 2d0*IB(n)*[conjg(ei2p(n))**2,i1]
  end do
  deallocate(ei2p)

  !* convolution
  allocate(Al(2,npix))
  do i = 1, 2
    call dft(A1(i,:),nn,D,-1)
    call dft(A2(i,:),nn,D,-1)
    Al(i,:) = A2(i,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
    Al(i,:) = dble(Al(i,:))
  end do
  deallocate(A1,A2)

  !* normalization
  do n = 1, npix
    if (present(eL).and.els(n)<eL(1).or.els(n)>eL(2))  cycle
    iAl = sum(Al(1:2,n))
    if(iAl<=0d0) cycle
    if (present(Ala))  Ala(n) = 1d0 / iAl
  end do

  deallocate(Al)

end subroutine alflat_eb_cb_nobb


end module rotrecflat


