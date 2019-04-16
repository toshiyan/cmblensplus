!///////////////////////////////////////////////////////////!
! * Tau Reconstruction Kernel
!///////////////////////////////////////////////////////////!

module taurecflat
  use myconst, only: iu, i1, pi, dlc, twopi
  use anaflat, only: elarrays
  use myfftw,  only: dft
  implicit none

  private iu, i1, pi, dlc, twopi
  private elarrays
  private dft

contains 


subroutine quadte_pr(nn,D,T,E,TE,rL,mlm)
  implicit none
!
! [input]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   TE(:)   --- TE correlation in 2D grids
!   T(:)    --- Temperature
!   E(:)    --- E-modes
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: TE(:), D(1:2)
  complex(dlc), intent(in) :: T(:), E(:)
!
! (optional)
  complex(dlc), intent(out), optional :: mlm(:)
!
! [internal]
  integer :: i, j, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: aT(:), aE(:), ei2p(:), bT(:,:), bE(:,:), alm(:,:)

  npix = nn(1)*nn(2)
  allocate(aE(npix),aT(npix),bE(3,npix),bT(3,npix),l(2,npix),els(npix),ei2p(npix)); aE = 0d0; bE = 0d0; aT = 0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,ei2p=ei2p)

  if (size(TE)/=npix)  stop 'error (quadte): size of TE /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    aE(n)     = E(n)*conjg(ei2p(n))
    bE(3,n)   = TE(n)*T(n)*ei2p(n)
    aT(n)     = T(n)
    bT(3,n)   = TE(n)*E(n)
  end do 

  !* convolution
  call dft(aE,nn,D,-1)
  call dft(aT,nn,D,-1)
  call dft(bE(3,:),nn,D,-1)
  call dft(bT(3,:),nn,D,-1)
  allocate(alm(3,npix))
  alm(3,:) = aimag(aE*bE(3,:)) - aT*bT(3,:)
  deallocate(aE,bE,aT,bT)
  call dft(alm(3,:),nn,D,1)

  !* form estimator 
  if(present(mlm)) mlm = alm(3,:)
  deallocate(alm,l,els)

end subroutine quadte_pr


subroutine quadtb_pr(nn,D,T,B,TE,rL,mlm)
  implicit none
!
! [input]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   TE(:)   --- TE correlation in 2D grids
!   T(:)    --- Temperature
!   B(:)    --- B-modes
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: TE(:), D(1:2)
  complex(dlc), intent(in) :: T(:), B(:)
!
! (optional)
  complex(dlc), intent(out), optional :: mlm(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:)
  complex(dlc), allocatable :: wB(:),wT(:),alm(:),ei2p(:)

  npix = nn(1)*nn(2)
  allocate(wB(npix),wT(npix),l(2,npix),els(npix),ei2p(npix)); wB=0d0; wT=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,ei2p=ei2p)

  if (size(TE)/=npix)  stop 'error (quadtb): size of TE /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wB(n) = B(n)*conjg(ei2p(n))
    wT(n) = TE(n)*T(n)*ei2p(n)
  end do 

  !* convolution
  call dft(wB,nn,D,-1)
  call dft(wT,nn,D,-1)
  allocate(alm(npix))
  alm = aimag(wB*wT)
  deallocate(wB,wT)
  call dft(alm,nn,D,1)

  !* form estimator 
  if(present(mlm)) mlm = alm
  deallocate(alm,l,els)

end subroutine quadtb_pr


subroutine quadeb_pr(nn,D,E,B,EE,rL,tlm)
!
  implicit none
! [inputs]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   EE(:)   --- E-mode Cl in 2D grids
!   E(:)    --- filtered E-modes
!   B(:)    --- filtered B-modes
!   eL(1:2) --- multipole range of output
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), D(1:2)
  complex(dlc), intent(in) :: E(:), B(:)
!
! [output]
!   tlm(:)  --- tau mode
  complex(dlc), intent(out) :: tlm(:)
!
! [internal]
  integer :: i, n, npix
  double precision :: els(nn(1)*nn(2))
  complex(dlc), allocatable :: wE(:), wB(:), alm(:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(wB(npix),wE(npix),ei2p(npix)); wB=0d0; wE=0d0
  call elarrays(nn,D,els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (quadeb): size(EE) /= npix'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wE(n) = EE(n)*E(n)*ei2p(n)
    wB(n) = B(n)*conjg(ei2p(n))
  end do
  deallocate(ei2p)

  !* convolution
  call dft(wB,nn,D,-1)
  call dft(wE,nn,D,-1)
  allocate(alm(npix))
  alm = aimag(wB*wE)
  deallocate(wB,wE)
  call dft(alm,nn,D,1)

  !* estimator 
  tlm = alm
  deallocate(alm)

end subroutine quadeb_pr


subroutine alflat_eb_pr(nn,D,fE,fB,CE,rL,eL,Alg)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CE(:), fE(:), fB(:), D(1:2)
  double precision, intent(out), optional :: Alg(:)
! [internal]
  integer :: n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc) :: vec(1:2)
  complex(dlc), allocatable :: Al(:,:), A1(:,:), A2(:,:), ei2p(:)

  write(*,*) 'EB flat pr'

  npix = nn(1)*nn(2)
  allocate(A1(3,npix),A2(3,npix),ei2p(npix));  A1=0d0;  A2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els,ei2p=ei2p)

  if (size(CE)/=npix)  stop 'error (alflat_eb_pr): CE is strange.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(1,n)  = CE(n)**2*fE(n)*ei2p(n)**2
    A1(2,n)  = CE(n)**2*fE(n)
    A2(1,n)  = - fB(n)*conjg(ei2p(n)**2)
    A2(2,n)  = fB(n)
  end do
  deallocate(ei2p)

  !* convolution
  allocate(Al(2,npix))
  do n = 1, 2
    call dft(A1(n,:),nn,D,-1)
    call dft(A2(n,:),nn,D,-1)
    Al(n,:) = A1(n,:)*A2(n,:)
    if (n==1)  Al(n,:) = dble(Al(n,:))
    call dft(Al(n,:),nn,D,1)
  end do
  deallocate(A1,A2)

  !* normalization
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    iAlg = sum(Al(1:2,n))*0.5d0
    if(iAlg<=0d0) cycle
    if (present(Alg))  Alg(n) = 1d0 / iAlg
  end do

  deallocate(Al)

end subroutine alflat_eb_pr


!/////////////////////////////////////////////////////////////////////////////////////!
! other quadratic estimators
!

subroutine quadtt_bhe(nn,D,T1,T2,TT,tL,mlm,slm)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), tL(1:2)
  double precision, intent(in) :: TT(:), D(1:2)
  complex(dlc), intent(in) :: T1(:), T2(:)
  complex(dlc), intent(out), optional :: mlm(:), slm(:)
  !internal
  integer :: i, n, npix
  double precision :: iC, els(nn(1)*nn(2))
  complex(dlc), allocatable :: aT(:), alm(:,:)

  call elarrays(nn,D,els=els)

  npix = nn(1)*nn(2)
  allocate(aT(npix),alm(2,npix)); aT=0d0; alm=0d0

  !* filtering
  do n = 1, npix
    if(tL(1)>els(n).or.els(n)>tL(2)) cycle
    aT(n)    = T1(n)
    alm(1,n) = TT(n)*T2(n)
    alm(2,n) = 0.5d0*T2(n)
  end do 

  !* convolution
  call dft(aT,nn,D,-1)
  do i = 1, 2
    call dft(alm(i,:),nn,D,-1)
    alm(i,:) = aT(:)*alm(i,:)
    call dft(alm(i,:),nn,D,1)
  end do
  deallocate(aT)

  !* form estimator
  if(present(mlm)) mlm = alm(1,:)
  if(present(slm)) slm = alm(2,:)

  deallocate(alm)

end subroutine quadtt_bhe


subroutine quadee_bhe(nn,D,E1,E2,EE,rL,mlm,slm)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rL(2)
  double precision, intent(in) :: EE(:), D(2)
  complex(dlc), intent(in) :: E1(:), E2(:)
  complex(dlc), intent(out), optional :: mlm(:),slm(:)
  !internal
  integer :: i, n, npix
  double precision :: els(nn(1)*nn(2))
  complex(dlc) :: ei2p(nn(1)*nn(2))
  complex(dlc), allocatable :: wE1(:), wE2(:,:), alm(:,:)

  call elarrays(nn,D,els=els,ei2p=ei2p)
  npix = nn(1)*nn(2)

  !* filtering
  allocate(wE1(npix),wE2(2,npix)); wE1=0d0; wE2=0d0
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wE1(n)   = E1(n)*conjg(ei2p(n))
    wE2(1,n) = EE(n)*E2(n)*ei2p(n)
    wE2(2,n) = E2(n)*ei2p(n)
  end do

  !* convolution
  call dft(wE1,nn,D,-1)
  call dft(wE2(1,:),nn,D,-1)
  call dft(wE2(2,:),nn,D,-1)
  allocate(alm(2,npix))
  alm(1,:) = real(wE1*wE2(1,:))
  alm(2,:) = real(wE1*wE2(2,:))
  deallocate(wE1,wE2)
  call dft(alm(1,:),nn,D,1)
  call dft(alm(2,:),nn,D,1)

  !* form estimator 
  if(present(mlm)) mlm = alm(1,:)
  if(present(slm)) slm = alm(2,:)
  deallocate(alm)

end subroutine quadee_bhe


subroutine quadeb_bhe(nn,D,E,B,EE,BB,rL,mlm,slm)
! * BHE for masking, point source, and polarization angle
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rL(2)
  double precision, intent(in) :: EE(:), BB(:), D(2)
  complex(dlc), intent(in) :: E(:),B(:)
  complex(dlc), intent(out), optional :: mlm(:),slm(:)
  !internal
  integer :: i, n, npix
  double precision :: els(nn(1)*nn(2))
  complex(dlc), allocatable :: X(:,:), alm(:), fE(:,:), fB(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(X(2,npix),fE(2,npix),fB(2,npix),ei2p(npix)); X=0d0; fE=0d0; fB=0d0
  call elarrays(nn,D,els=els,ei2p=ei2p)

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    X(:,n)  = [E(n),B(n)]*conjg(ei2p(n))
    fE(1,n) = EE(n)*E(n)*ei2p(n)
    fB(1,n) = BB(n)*B(n)*ei2p(n)
    fE(2,n) = E(n)*ei2p(n)*0.5d0
    fB(2,n) = B(n)*ei2p(n)*0.5d0
  end do 

  !* convolution
  do i = 1, 2
    call dft(X(i,:),nn,D,-1)
    call dft(fE(i,:),nn,D,-1)
    call dft(fB(i,:),nn,D,-1)
  end do

  !* estimator 
  allocate(alm(npix));  alm=0d0
  if(present(mlm)) then
    alm = aimag(X(2,:)*fE(1,:)-X(1,:)*fB(1,:))
    call dft(alm,nn,D,1)
    mlm = alm
  end if
  if(present(slm)) then
    alm = aimag(X(2,:)*fE(2,:)-X(1,:)*fB(2,:))
    call dft(alm,nn,D,1)
    slm = alm
  end if
  deallocate(X,fE,fB,alm)

end subroutine quadeb_bhe


end module taurecflat


