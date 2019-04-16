!///////////////////////////////////////////////////////////!
! * Lensing Reconstruction Kernel
!///////////////////////////////////////////////////////////!

module recflat
  use constants, only: iu, i1, pi, dlc, twopi
  use anaflat, only: elarrays, elxy
  use fftw,  only: dft
  implicit none

  interface quadtt
    module procedure quadtt_sym, quadtt_asym
  end interface quadtt

  interface alflat_tt
    module procedure alflat_tt_1d, alflat_tt_2d
  end interface alflat_tt

  private iu, i1, pi, dlc, twopi
  private elarrays, elxy
  private dft

contains 


subroutine quadtt_sym(nn,D,T,TT,rL,glm,clm,gtype)
  implicit none
!
! [input]
!   nn(1:2)    --- 2D grids
!   D(1:2)     --- map size (radian)
!   rL(1:2)    --- multipole range of reconstruction
!   TT(1:npix) --- TT correlation in 2D grids
!   T(1:npix)  --- Temperature
  integer, intent(in) :: nn(1:2)
  integer, intent(in) :: rL(1:2)
  double precision, intent(in) :: TT(:), D(1:2)
  complex(dlc), intent(in) :: T(:)
!
! (optional)
  character(*), intent(in), optional :: gtype
  complex(dlc), intent(out), optional :: glm(:),clm(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: aT(:), alm(:,:)

  npix = nn(1)*nn(2)
  allocate(aT(npix),alm(2,npix),l(2,npix),li(npix),els(npix));  aT=0d0;  alm=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),eli=li,els=els)

  !* kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  if (size(TT)/=npix)  stop 'error (quadtt): size of TT /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    aT(n)    = T(n)
    alm(:,n) = l(:,n)*TT(n)*T(n)
  end do

  !* convolution
  call dft(aT,nn,D,-1)
  do i = 1, 2
    call dft(alm(i,:),nn,D,-1)
    alm(i,:) = aT*alm(i,:)
    call dft(alm(i,:),nn,D,1)
  end do

  deallocate(aT)

  !* estimator 
  if(present(glm)) glm = sum(alm*l,dim=1)*li
  if(present(clm)) clm = (alm(1,:)*l(2,:)-alm(2,:)*l(1,:))*li
  deallocate(alm,els,l,li)

end subroutine quadtt_sym


subroutine quadtt_asym(nn,D,T1,T2,TT,rL,glm,clm,gtype)
  implicit none
! [input]
!   nn(1:2) --- 2D grids
!   D(1:2)  --- map size (radian)
!   rL(1:2) --- multipole range of reconstruction
!   TT(:)   --- TT correlation in 2D grids
!   T1(:)   --- Temperature 1
!   T2(:)   --- Temperature 2
  integer, intent(in) :: nn(1:2)
  integer, intent(in) :: rL(1:2)
  double precision, intent(in) :: TT(:), D(1:2)
  complex(dlc), intent(in) :: T1(:), T2(:)
!
! (optional)
  character(*), intent(in), optional :: gtype
  complex(dlc), intent(out), optional :: glm(:),clm(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: aT(:), alm(:,:), blm(:,:)

  npix = nn(1)*nn(2)
  allocate(aT(npix),alm(2,npix),blm(2,npix),l(2,npix),li(npix),els(npix));  aT=0d0;  alm=0d0;  blm=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),eli=li,els=els)

  !* kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  if (size(TT)/=npix)  stop 'error (quadtt): size of TT /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    aT(n)    = T1(n)
    alm(:,n) = 0.5d0*l(:,n)*TT(n)*T2(n)
  end do

  !* convolution
  call dft(aT,nn,D,-1)
  do i = 1, 2
    call dft(alm(i,:),nn,D,-1)
    alm(i,:) = aT*alm(i,:)
    call dft(alm(i,:),nn,D,1)
  end do

  !* filtering
  aT = 0d0
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    aT(n)    = T2(n)
    blm(:,n) = 0.5d0*l(:,n)*TT(n)*T1(n)
  end do

  !* convolution
  call dft(aT,nn,D,-1)
  do i = 1, 2
    call dft(blm(i,:),nn,D,-1)
    blm(i,:) = aT*blm(i,:)
    call dft(blm(i,:),nn,D,1)
  end do

  deallocate(aT)

  !* estimator 
  if(present(glm)) glm = (sum(alm*l,dim=1)+sum(blm*l,dim=1))*li
  if(present(clm)) clm = (alm(1,:)*l(2,:)-alm(2,:)*l(1,:)+blm(1,:)*l(2,:)-blm(2,:)*l(1,:))*li
  deallocate(alm,blm,els,l,li)

end subroutine quadtt_asym


subroutine quadte(nn,D,T,E,TE,rL,glm,clm,gtype)
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
  character(*), intent(in), optional :: gtype
  complex(dlc), intent(out), optional :: glm(:),clm(:)
!
! [internal]
  integer :: i, j, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: aT(:,:), bT(:), aE(:), bE(:,:), alm(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(aT(2,npix),bT(npix),aE(npix),bE(2,npix),l(2,npix),li(npix),els(npix),ei2p(npix)); aE=0d0; bE=0d0; aT=0d0; bT=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),eli=li,els=els,ei2p=ei2p)

  !* kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  if (size(TE)/=npix)  stop 'error (quadte): size of TE /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    aT(1:2,n) = l(1:2,n)*TE(n)*T(n)*ei2p(n)
    aE(n)     = E(n)*conjg(ei2p(n))
    bT(n)     = T(n)
    bE(1:2,n) = l(1:2,n)*TE(n)*E(n)
  end do 

  !* convolution
  call dft(aE,nn,D,-1)
  call dft(bT,nn,D,-1)
  do i = 1, 2
    call dft(aT(i,:),nn,D,-1)
    call dft(bE(i,:),nn,D,-1)
  end do
  allocate(alm(2,npix)); alm=0d0
  alm(1,:) = iu*aimag(aE*aT(1,:)) + bT*bE(1,:)
  alm(2,:) = iu*aimag(aE*aT(2,:)) + bT*bE(2,:)
  !alm(1,:) = bT*bE(1,:)
  !alm(2,:) = bT*bE(2,:)
  deallocate(aE,bE,aT,bT)
  do i = 1, 2
    call dft(alm(i,:),nn,D,1)
  end do

  !* form estimator 
  if(present(glm)) glm = (alm(1,:)*l(1,:)+alm(2,:)*l(2,:))*li
  if(present(clm)) clm = (alm(1,:)*l(2,:)-alm(2,:)*l(1,:))*li
  deallocate(alm,l,els,li)

end subroutine quadte


subroutine quadtb(nn,D,T,B,TE,rL,glm,clm,gtype)
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
  character(*), intent(in), optional :: gtype
  complex(dlc), intent(out), optional :: glm(:),clm(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: wB(:),wT(:,:),alm(:,:),ei2p(:)

  npix = nn(1)*nn(2)
  allocate(wB(npix),wT(2,npix),l(2,npix),li(npix),els(npix),ei2p(npix)); wB=0d0; wT=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),eli=li,els=els,ei2p=ei2p)

  !* kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  if (size(TE)/=npix)  stop 'error (quadtb): size of TE /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wB(n)     = B(n)*conjg(ei2p(n))
    wT(1:2,n) = l(1:2,n)*TE(n)*T(n)*ei2p(n)
  end do 

  !* convolution
  call dft(wB,nn,D,-1)
  do i = 1, 2
    call dft(wT(i,:),nn,D,-1)
  end do
  allocate(alm(2,npix))
  alm(1,:) = real(wB(:)*wT(1,:))/iu
  alm(2,:) = real(wB(:)*wT(2,:))/iu
  deallocate(wB,wT)
  do i = 1, 2
    call dft(alm(i,:),nn,D,1)
  end do

  !* form estimator 
  if(present(glm)) glm = (alm(1,:)*l(1,:)+alm(2,:)*l(2,:))*li
  if(present(clm)) clm = (alm(1,:)*l(2,:)-alm(2,:)*l(1,:))*li
  deallocate(alm,l,li,els)

end subroutine quadtb


subroutine quadee(nn,D,E1,E2,EE,rL,glm,clm,eL,gtype)
  implicit none
! [input]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   fC(:)   --- E-mode Cl in 2D grids
!   E1(:)   --- E-modes 1
!   E2(:)   --- E-modes 2
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), D(1:2)
  complex(dlc), intent(in) :: E1(:), E2(:)
!
! (optional)
!   eL(1:2) --- multipole range of output
!   glm(:)  --- gradient mode
!   clm(:)  --- curl mode
  character(*), intent(in), optional :: gtype
  integer, intent(in), optional :: eL(1:2)
  complex(dlc), intent(out), optional :: glm(:), clm(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: wE1(:), wE2(:,:), alm(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(wE1(npix),wE2(2,npix),l(2,npix),li(npix),els(npix),ei2p(npix)); wE1=0d0; wE2=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),eli=li,els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (quadee): size of EE /= npix.'

  !* kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wE1(n)   = E1(n)*conjg(ei2p(n))
    wE2(:,n) = l(:,n)*EE(n)*E2(n)*ei2p(n)
  end do 

  !* convolution
  call dft(wE1,nn,D,-1)
  call dft(wE2(1,:),nn,D,-1)
  call dft(wE2(2,:),nn,D,-1)
  allocate(alm(2,npix))
  alm(1,:) = aimag(wE1*wE2(1,:))*iu
  alm(2,:) = aimag(wE1*wE2(2,:))*iu
  deallocate(wE1,wE2)
  call dft(alm(1,:),nn,D,1)
  call dft(alm(2,:),nn,D,1)

  !* estimator 
  do n = 1, npix
    if(present(eL).and.(els(n)<eL(1).or.els(n)>eL(2))) alm(:,n)=0d0
    if(present(glm)) glm(n) = (l(1,n)*alm(1,n)+l(2,n)*alm(2,n))*li(n)
    if(present(clm)) clm(n) = (l(2,n)*alm(1,n)-l(1,n)*alm(2,n))*li(n)
  end do
  deallocate(alm,l,els,li)

end subroutine quadee


subroutine quadeb(nn,D,E,B,EE,rL,glm,clm,eL,gtype)
  implicit none
! [inputs]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   fC(:)   --- E-mode Cl in 2D grids
!   E(:)    --- filtered E-modes
!   B(:)    --- filtered B-modes
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), D(1:2)
  complex(dlc), intent(in) :: E(:), B(:)
!
! (optional)
!   gtype   --- phi or kappa
!   eL(1:2) --- multipole range of output
!   glm(:)  --- gradient mode
!   clm(:)  --- curl mode
  character(*), intent(in), optional :: gtype
  integer, intent(in), optional :: eL(1:2)
  complex(dlc), intent(out), optional :: glm(:), clm(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: wB(:), alm(:,:), wE(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(wB(npix),wE(2,npix),ei2p(npix),l(2,npix),els(npix),li(npix));  wB=0d0;  wE=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,eli=li,ei2p=ei2p)

  !* kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  if (size(EE)/=npix)  stop 'error (quadeb): size of EE /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wB(n)   = B(n)*conjg(ei2p(n))
    wE(:,n) = l(:,n)*EE(n)*E(n)*ei2p(n)
  end do
  deallocate(ei2p)

  !* convolution
  call dft(wB,nn,D,-1)
  call dft(wE(1,:),nn,D,-1)
  call dft(wE(2,:),nn,D,-1)
  allocate(alm(2,npix))
  alm(1,:) = dble(wB*wE(1,:))
  alm(2,:) = dble(wB*wE(2,:))
  deallocate(wB,wE)
  call dft(alm(1,:),nn,D,1)
  call dft(alm(2,:),nn,D,1)

  !* estimator 
  do n = 1, npix
    if(present(eL).and.(els(n)<eL(1).or.els(n)>eL(2))) alm(:,n)=0d0
    if(present(glm)) glm(n) = sum(l(:,n)*alm(:,n))/iu * li(n)
    if(present(clm)) clm(n) = (l(2,n)*alm(1,n)-l(1,n)*alm(2,n))/iu * li(n)
  end do

  deallocate(alm,l,els,li)

end subroutine quadeb


subroutine quadbb(nn,D,B1,B2,BB,rL,glm,clm,mlm,slm)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rL(2)
  double precision, intent(in) :: BB(:), D(2)
  complex(dlc), intent(in) :: B1(:), B2(:)
  complex(dlc), intent(out), optional :: glm(:),clm(:),mlm(:),slm(:)
  !internal
  integer :: i, n, npix
  double precision :: l(2,nn(1)*nn(2)), els(nn(1)*nn(2))
  complex(dlc) ::  ei2p(nn(1)*nn(2))
  complex(dlc), allocatable :: wB1(:), wB2(:,:), alm(:,:)

  npix = nn(1)*nn(2)
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,ei2p=ei2p)

  if (size(BB)/=npix)  stop 'error (quadbb): size of BB is not npix.'

  allocate(wB1(npix),wB2(3,npix)); wB1=0d0; wB2=0d0
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wB1(n)     = B1(n)*conjg(ei2p(n))
    wB2(3,n)   = BB(n)*B2(n)*ei2p(n)
    wB2(1:2,n) = l(1:2,n)*wB2(3,n)
  end do 

  !* convolution
  call dft(wB1,nn,D,-1)
  do i = 1, 4
    call dft(wB2(i,:),nn,D,-1)
  end do
  allocate(alm(4,npix))
  alm(1,:) = aimag(wB1*wB2(1,:))*iu
  alm(2,:) = aimag(wB1*wB2(2,:))*iu
  alm(3,:) = real(wB1*wB2(3,:))
  deallocate(wB1,wB2)
  do i = 1, 4
    call dft(alm(i,:),nn,D,1)
  end do

  !* form estimator 
  if(present(glm)) glm = -(alm(1,:)*l(1,:)+alm(2,:)*l(2,:))
  if(present(clm)) clm = -(alm(1,:)*l(2,:)-alm(2,:)*l(1,:))
  if(present(mlm)) mlm = alm(3,:)
  if(present(slm)) slm = alm(4,:)
  deallocate(alm)

end subroutine quadbb


!/////////////////////////////////////////////////////////////////////////////////////!
! 2D normalization

subroutine alflat_tt_1d(nn,D,OC,CT,rL,eL,Alg,Alc)
  implicit none

  !I/O
  ! OC(nn(1)*nn(2)) --- inverse of observed cl in 2D
  ! CT(nn(1)*nn(2)) --- lensed cl in 2D
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CT(:), OC(:), D(1:2)

  !(optional)
  ! Alg(nn(1)*nn(2)) --- normalization (gradient)
  ! Alc(nn(1)*nn(2)) --- normalization (curl)
  double precision, intent(out), optional :: Alg(:), Alc(:)

  !internal
  integer :: i, n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc) :: vec
  complex(dlc), allocatable :: Al(:,:), Bl(:,:), A1(:,:), A2(:,:), B1(:,:), B2(:,:)

  npix = nn(1)*nn(2)
  allocate(A1(3,npix),A2(1,npix),B1(2,npix),B2(2,npix)); A1=0d0; A2=0d0; B1=0d0; B2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els)

  if (size(CT)/=npix)  stop 'error (alflat_tt): size of TT is not npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    vec     = CT(n)**2 * OC(n)
    A1(1,n) = ll(1,n)**2 * vec
    A1(2,n) = 2d0*ll(1,n)*ll(2,n) * vec
    A1(3,n) = ll(2,n)**2 * vec
    A2(1,n) = OC(n)
    B1(:,n) = ll(:,n) * CT(n) * OC(n)
    B2(:,n) = B1(:,n)
  end do

  !* convolution
  call dft(A2(1,:),nn,D,-1)
  allocate(Al(3,npix))
  do i = 1, 3
    call dft(A1(i,:),nn,D,-1)
    Al(i,:) = A2(1,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
  end do
  deallocate(A1,A2)

  do i = 1, 2
    call dft(B1(i,:),nn,D,-1)
    call dft(B2(i,:),nn,D,-1)
  end do
  allocate(Bl(4,npix))
  Bl(1,:) = B1(1,:)*B2(1,:)
  Bl(2,:) = B1(2,:)*B2(2,:)
  Bl(3,:) = B1(1,:)*B2(2,:)
  Bl(4,:) = B1(2,:)*B2(1,:)
  do i = 1, 4
    call dft(Bl(i,:),nn,D,1)
  end do
  deallocate(B1,B2)

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    if (present(Alg)) then
      iAlg = ll(1,n)**2*Al(1,n) + ll(1,n)*ll(2,n)*Al(2,n) + ll(2,n)**2*Al(3,n) &
        + ll(1,n)**2*Bl(1,n) + ll(2,n)**2*Bl(2,n) + ll(1,n)*ll(2,n)*sum(Bl(3:4,n))
      if(iAlg>0d0) Alg(n) = 1d0/iAlg
    end if
    if (present(Alc)) then
      iAlc = ll(2,n)**2*Al(1,n) - ll(1,n)*ll(2,n)*Al(2,n) + ll(1,n)**2*Al(3,n) &
        + ll(2,n)**2*Bl(1,n) + ll(1,n)**2*Bl(2,n) - ll(1,n)*ll(2,n)*sum(Bl(3:4,n))
      if(iAlc>0d0) Alc(n) = 1d0/iAlc
    end if
  end do

  deallocate(Al)

end subroutine alflat_tt_1d


subroutine alflat_tt_2d(nn,D,OC,CT,rL,eL,Alg,Alc)
  implicit none

  !I/O
  ! OC(nn(1),nn(2)) --- inverse of observed cl in 2D
  ! CT(nn(1),nn(2)) --- lensed cl in 2D
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CT(:,:), OC(:,:), D(1:2)

  !(optional)
  ! Alg(nn(1),nn(2)) --- normalization (gradient)
  ! Alc(nn(1),nn(2)) --- normalization (curl)
  double precision, intent(out), optional :: Alg(:,:), Alc(:,:)

  !internal
  integer :: i, j
  double precision :: fac, lx, ly, ll, iAlg, iAlc
  complex(dlc), allocatable :: Al(:,:,:), Bl(:,:,:), A1(:,:,:), A2(:,:,:), B(:,:,:)

  allocate(A1(3,nn(1),nn(2)),A2(1,nn(1),nn(2)),B(2,nn(1),nn(2))); A1=0d0; A2=0d0; B=0d0

  !* filtering
  do i = 1, nn(1)
    lx = elxy(i,nn(1),D(1))
    do j = 1, nn(2)
      ly = elxy(j,nn(2),D(2))
      ll = dsqrt(lx**2+ly**2)
      if(rL(1)>ll.or.ll>rL(2)) cycle
      fac       = CT(i,j)**2 * OC(i,j)
      A1(1,i,j) = lx**2 * fac
      A1(2,i,j) = lx*ly * fac
      A1(3,i,j) = ly**2 * fac
      A2(1,i,j) = OC(i,j)
      B(1,i,j)  = lx * CT(i,j) * OC(i,j)
      B(2,i,j)  = ly * CT(i,j) * OC(i,j)
    end do
  end do

  !* convolution
  call dft(A2(1,:,:),nn,D,-1)
  allocate(Al(3,nn(1),nn(2)))
  do i = 1, 3
    call dft(A1(i,:,:),nn,D,-1)
    Al(i,:,:) = A2(1,:,:)*A1(i,:,:)
    call dft(Al(i,:,:),nn,D,1)
  end do
  deallocate(A1,A2)

  do i = 1, 2
    call dft(B(i,:,:),nn,D,-1)
  end do
  allocate(Bl(3,nn(1),nn(2)))
  Bl(1,:,:) = B(1,:,:)**2
  Bl(2,:,:) = B(2,:,:)**2
  Bl(3,:,:) = B(1,:,:)*B(2,:,:)
  do i = 1, 3
    call dft(Bl(i,:,:),nn,D,1)
  end do
  deallocate(B)

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0

  do i = 1, nn(1)
    lx = elxy(i,nn(1),D(1))
    do j = 1, nn(2)
      ly = elxy(j,nn(2),D(2))
      ll = dsqrt(lx**2+ly**2)
      if (ll<eL(1).or.ll>eL(2))  cycle
      if (present(Alg)) then
        iAlg = lx**2*Al(1,i,j) + 2d0*lx*ly*Al(2,i,j) + ly**2*Al(3,i,j) + lx**2*Bl(1,i,j) + ly**2*Bl(2,i,j) + 2d0*lx*ly*Bl(3,i,j)
        if(iAlg>0d0) Alg(i,j) = 1d0/iAlg
      end if
      if (present(Alc)) then
        iAlc = ly**2*Al(1,i,j) - 2d0*lx*ly*Al(2,i,j) + lx**2*Al(3,i,j) + ly**2*Bl(1,i,j) + lx**2*Bl(2,i,j) - 2d0*lx*ly*Bl(3,i,j)
        if(iAlc>0d0) Alc(i,j) = 1d0/iAlc
      end if
    end do
  end do

  deallocate(Al,Bl)

end subroutine alflat_tt_2d


subroutine alflat_te(nn,D,fT,fE,TE,rL,eL,Alg,Alc)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), eL(1:2), rL(1:2)
  double precision, intent(in) :: TE(:), fT(:), fE(:), D(1:2)
  double precision, intent(out), optional :: Alg(:), Alc(:)
! [internal]
  integer :: i, n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc), allocatable :: ei2p(:), C1(:,:), C2(:,:), C3(:,:), A1(:,:), A2(:,:), A3(:), B1(:,:), B2(:,:), B3(:,:)

  npix = nn(1)*nn(2)
  allocate(A1(6,npix),A2(2,npix),A3(npix),B1(2,npix),B2(2,npix),B3(3,npix),ei2p(npix)); A1=0d0; A2=0d0; A3=0d0; B1=0d0; B2=0d0; B3=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els,ei2p=ei2p)

  if (size(TE)/=npix)  stop 'error (alflat_te): size of TE is not npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(1,n) = fT(n) * TE(n)**2 * ei2p(n)**2/2d0 * ll(1,n)**2
    A1(2,n) = fT(n) * TE(n)**2 * ei2p(n)**2/2d0 * ll(1,n)*ll(2,n)*2d0
    A1(3,n) = fT(n) * TE(n)**2 * ei2p(n)**2/2d0 * ll(2,n)**2
    A1(4,n) = fT(n) * TE(n)**2 / 2d0 * ll(1,n)**2
    A1(5,n) = fT(n) * TE(n)**2 / 2d0 * ll(1,n)*ll(2,n)*2d0
    A1(6,n) = fT(n) * TE(n)**2 / 2d0 * ll(2,n)**2
    A2(1,n) = fT(n) * TE(n) * ei2p(n) * ll(1,n)
    A2(2,n) = fT(n) * TE(n) * ei2p(n) * ll(2,n)
    A3(n)   = fT(n)
    B1(1,n) = fE(n) * conjg(ei2p(n))**2
    B1(2,n) = fE(n)
    B2(1,n) = fE(n) * TE(n) * conjg(ei2p(n)) * ll(1,n)
    B2(2,n) = fE(n) * TE(n) * conjg(ei2p(n)) * ll(2,n)
    B3(1,n) = fE(n) * TE(n)**2 * ll(1,n)**2
    B3(2,n) = fE(n) * TE(n)**2 * ll(1,n)*ll(2,n)*2d0
    B3(3,n) = fE(n) * TE(n)**2 * ll(2,n)**2
  end do
  deallocate(ei2p)

  !* inverse FT
  do i = 1, 6
    call dft(A1(i,:),nn,D,-1)
  end do
  do i = 1, 3
    call dft(B3(i,:),nn,D,-1)
  end do
  do i = 1, 2
    call dft(A2(i,:),nn,D,-1)
    call dft(B1(i,:),nn,D,-1)
    call dft(B2(i,:),nn,D,-1)
  end do
  call dft(A3,nn,D,-1)

  !* convolution
  allocate(C1(3,npix),C2(3,npix),C3(3,npix))
  do i = 1, 3
    C1(i,:) = A1(i,:)*B1(1,:) + A1(i+3,:)*B1(2,:)
    call dft(C1(i,:),nn,D,1)
    C1(i,:) = dble(C1(i,:))
  end do
  C2(1,:) = A2(1,:)*B2(1,:)
  C2(2,:) = A2(1,:)*B2(2,:) + A2(2,:)*B2(1,:)
  C2(3,:) = A2(2,:)*B2(2,:)
  do i = 1, 3
    call dft(C2(i,:),nn,D,1)
    C2(i,:) = dble(C2(i,:))
  end do
  do i = 1, 3
    C3(i,:) = A3*B3(i,:)
    call dft(C3(i,:),nn,D,1)
  end do
  deallocate(A1,A2,A3,B1,B2,B3)

  C1 = C1 + 2*C2 + C3

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    if (present(Alg)) then
      iAlg = ll(1,n)**2*C1(1,n) + ll(2,n)**2*C1(3,n) + ll(1,n)*ll(2,n)*C1(2,n)
      if(iAlg>0d0) Alg(n) = 1d0/iAlg
    end if
    if (present(Alc)) then
      iAlc = ll(2,n)**2*C1(1,n) + ll(1,n)**2*C1(3,n) - ll(1,n)*ll(2,n)*C1(2,n)
      if(iAlc>0d0) Alc(n) = 1d0/iAlc
    end if
  end do
  deallocate(C1,C2,C3)

end subroutine alflat_te


subroutine alflat_tb(nn,D,fT,fB,CTE,rL,eL,Alg,Alc)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CTE(:), fT(:), fB(:), D(1:2)
  double precision, intent(out), optional :: Alg(:), Alc(:)
! [internal]
  integer :: i, n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc) :: vec(1:2)
  complex(dlc), allocatable :: Al(:,:), A1(:,:), A2(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(A1(6,npix),A2(2,npix),ei2p(npix));  A1=0d0;  A2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els,ei2p=ei2p)

  if (size(CTE)/=npix)  stop 'error (alflat_tb): CTE is strange.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    vec       = [-ei2p(n)**2,i1] * CTE(n)**2 * fT(n)
    A1(1:2,n) = ll(1,n)**2 * vec
    A1(3:4,n) = 2d0*ll(1,n)*ll(2,n) * vec
    A1(5:6,n) = ll(2,n)**2 * vec
    A2(1:2,n) = 0.5d0*fB(n)*[conjg(ei2p(n))**2,i1]
  end do
  deallocate(ei2p)

  !* convolution
  call dft(A2(1,:),nn,D,-1)
  call dft(A2(2,:),nn,D,-1)
  allocate(Al(6,npix))
  do i = 1, 6
    call dft(A1(i,:),nn,D,-1)
    if (mod(i,2)==0)  Al(i,:) = A2(2,:)*A1(i,:)
    if (mod(i,2)==1)  Al(i,:) = A2(1,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
    if (mod(i,2)==1)  Al(i,:) = dble(Al(i,:))
  end do
  deallocate(A1,A2)

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    iAlg = ll(1,n)**2*sum(Al(1:2,n)) + ll(1,n)*ll(2,n)*sum(Al(3:4,n)) + ll(2,n)**2*sum(Al(5:6,n))
    iAlc = ll(2,n)**2*sum(Al(1:2,n)) - ll(1,n)*ll(2,n)*sum(Al(3:4,n)) + ll(1,n)**2*sum(Al(5:6,n))
    if(iAlg<=0d0.or.iAlc<=0d0) cycle
    if (present(Alg))  Alg(n) = 1d0 / iAlg
    if (present(Alc))  Alc(n) = 1d0 / iAlc
  end do

  deallocate(Al)

end subroutine alflat_tb


subroutine alflat_ee(nn,D,fE,EE,rL,eL,Alg,Alc)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), fE(:), D(1:2)
  double precision, intent(out), optional :: Alg(:), Alc(:)
! [internal]
  integer :: i, n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc) :: vec(1:2)
  complex(dlc), allocatable :: Al(:,:), Bl(:,:), A1(:,:), A2(:,:), ei2p(:), B1(:,:), B2(:,:)

  npix = nn(1)*nn(2)
  allocate(A1(6,npix),A2(2,npix),B1(4,npix),B2(4,npix),ei2p(npix)); A1=0d0; A2=0d0; B1=0d0; B2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (alflat_ee): size of EE is not npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    vec       = [ei2p(n)**2,i1] * EE(n)**2 * fE(n)
    A1(1:2,n) = ll(1,n)**2 * vec
    A1(3:4,n) = 2d0*ll(1,n)*ll(2,n) * vec
    A1(5:6,n) = ll(2,n)**2 * vec
    A2(1:2,n) = 0.5d0*fE(n)*[conjg(ei2p(n))**2,i1]
    vec       = [ei2p(n)**2,i1] * EE(n) * fE(n)
    B1(1:2,n) = ll(1,n) * vec
    B1(3:4,n) = ll(2,n) * vec
    B2(:,n)   = conjg(B1(:,n))*0.5d0
  end do
  deallocate(ei2p)

  !* convolution
  call dft(A2(1,:),nn,D,-1)
  call dft(A2(2,:),nn,D,-1)
  allocate(Al(6,npix))
  do i = 1, 6
    call dft(A1(i,:),nn,D,-1)
    if (mod(i,2)==0)  Al(i,:) = A2(2,:)*A1(i,:)
    if (mod(i,2)==1)  Al(i,:) = A2(1,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
    if (mod(i,2)==1)  Al(i,:) = dble(Al(i,:))
  end do
  deallocate(A1,A2)

  do i = 1, 4
    call dft(B1(i,:),nn,D,-1)
    call dft(B2(i,:),nn,D,-1)
  end do
  allocate(Bl(8,npix))
  do i = 1, 4
    Bl(i,:) = B1(i,:)*B2(i,:)
    if (i==1)  Bl(i+4,:) = B1(i,:)*B2(3,:)
    if (i==2)  Bl(i+4,:) = B1(i,:)*B2(4,:)
    if (i==3)  Bl(i+4,:) = B1(i,:)*B2(1,:)
    if (i==4)  Bl(i+4,:) = B1(i,:)*B2(2,:)
    call dft(Bl(i,:),nn,D,1)
    call dft(Bl(i+4,:),nn,D,1)
    if (mod(i,2)==1)  Bl(i,:) = dble(Bl(i,:))
    if (mod(i,2)==1)  Bl(i+4,:) = dble(Bl(i+4,:))
  end do
  deallocate(B1,B2)

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    if (present(Alg)) then
      iAlg = ll(1,n)**2*sum(Al(1:2,n)) + ll(1,n)*ll(2,n)*sum(Al(3:4,n)) + ll(2,n)**2*sum(Al(5:6,n)) &
        + ll(1,n)**2*sum(Bl(1:2,n)) + ll(2,n)**2*sum(Bl(3:4,n)) + ll(1,n)*ll(2,n)*sum(Bl(5:8,n))
      if(iAlg>0d0) Alg(n) = 1d0/iAlg
    end if
    if (present(Alc)) then
      iAlc = ll(2,n)**2*sum(Al(1:2,n)) - ll(1,n)*ll(2,n)*sum(Al(3:4,n)) + ll(1,n)**2*sum(Al(5:6,n)) &
        + ll(2,n)**2*sum(Bl(1:2,n)) + ll(1,n)**2*sum(Bl(3:4,n)) - ll(1,n)*ll(2,n)*sum(Bl(5:8,n))
      if(iAlc>0d0) Alc(n) = 1d0/iAlc
    end if
  end do

  deallocate(Al)

end subroutine alflat_ee


subroutine alflat_eb(nn,D,fE,fB,CE,rL,eL,Alg,Alc)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CE(:), fE(:), fB(:), D(1:2)
  double precision, intent(out), optional :: Alg(:), Alc(:)
! [internal]
  integer :: i, n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc) :: vec(1:2)
  complex(dlc), allocatable :: Al(:,:), A1(:,:), A2(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(A1(6,npix),A2(2,npix),ei2p(npix));  A1=0d0;  A2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els,ei2p=ei2p)

  if (size(CE)/=npix)  stop 'error (alflat_eb): CE is strange.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    vec       = [-ei2p(n)**2,i1] * CE(n)**2 * fE(n)
    A1(1:2,n) = ll(1,n)**2 * vec
    A1(3:4,n) = 2d0*ll(1,n)*ll(2,n) * vec
    A1(5:6,n) = ll(2,n)**2 * vec
    A2(1:2,n) = 0.5d0*fB(n)*[conjg(ei2p(n))**2,i1]
  end do
  deallocate(ei2p)

  !* convolution
  call dft(A2(1,:),nn,D,-1)
  call dft(A2(2,:),nn,D,-1)
  allocate(Al(6,npix))
  do i = 1, 6
    call dft(A1(i,:),nn,D,-1)
    if (mod(i,2)==0)  Al(i,:) = A2(2,:)*A1(i,:)
    if (mod(i,2)==1)  Al(i,:) = A2(1,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
    if (mod(i,2)==1)  Al(i,:) = dble(Al(i,:))
  end do
  deallocate(A1,A2)

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    iAlg = ll(1,n)**2*sum(Al(1:2,n)) + ll(1,n)*ll(2,n)*sum(Al(3:4,n)) + ll(2,n)**2*sum(Al(5:6,n))
    iAlc = ll(2,n)**2*sum(Al(1:2,n)) - ll(1,n)*ll(2,n)*sum(Al(3:4,n)) + ll(1,n)**2*sum(Al(5:6,n))
    if(iAlg<=0d0.or.iAlc<=0d0) cycle
    if (present(Alg))  Alg(n) = 1d0 / iAlg
    if (present(Alc))  Alc(n) = 1d0 / iAlc
  end do

  deallocate(Al)

end subroutine alflat_eb


! asymmetric case
subroutine dn0flat_tt(nn,D,OC1,OC2,CT,rL,eL,Alg,Alc)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CT(:), OC1(:), OC2(:), D(1:2)
  !optional
  double precision, intent(out), optional :: Alg(:), Alc(:)
  !internal
  integer :: i, n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc), allocatable :: Al(:,:), Bl(:,:), A1(:,:), A2(:,:), B1(:,:), B2(:,:)

  npix = nn(1)*nn(2)
  allocate(A1(4,npix),A2(4,npix),B1(2,npix),B2(2,npix)); A1=0d0; A2=0d0; B1=0d0; B2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els)

  if (size(CT)/=npix)  stop 'error (alflat_tt): size of TT is not npix.'

  !* for convolution
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(1,n) = OC1(n) * CT(n)**2 * ll(1,n)**2
    A1(2,n) = OC1(n) * CT(n)**2 * 2d0*ll(1,n)*ll(2,n)
    A1(3,n) = OC1(n) * CT(n)**2 * ll(2,n)**2
    A1(4,n) = OC2(n)
    A2(1,n) = OC2(n) * CT(n)**2 * ll(1,n)**2
    A2(2,n) = OC2(n) * CT(n)**2 * 2d0*ll(1,n)*ll(2,n)
    A2(3,n) = OC2(n) * CT(n)**2 * ll(2,n)**2
    A2(4,n) = OC1(n)
    B1(1:2,n) = ll(1:2,n) * CT(n) * OC1(n)
    B2(1:2,n) = ll(1:2,n) * CT(n) * OC2(n)
  end do

  !* convolution
  call dft(A1(4,:),nn,D,-1)
  call dft(A2(4,:),nn,D,-1)
  allocate(Al(6,npix))
  do i = 1, 3
    call dft(A1(i,:),nn,D,-1)
    Al(i,:) = A1(4,:)*A1(i,:)
    call dft(Al(i,:),nn,D,1)
    call dft(A2(i,:),nn,D,-1)
    Al(i+3,:) = A2(4,:)*A2(i,:)
    call dft(Al(i+3,:),nn,D,1)
  end do
  deallocate(A1,A2)

  do i = 1, 2
    call dft(B1(i,:),nn,D,-1)
    call dft(B2(i,:),nn,D,-1)
  end do
  allocate(Bl(4,npix))
  Bl(1:2,:) = B1(1:2,:)*B2(1:2,:)
  Bl(3,:)   = B1(1,:)*B2(2,:)
  Bl(4,:)   = B1(2,:)*B2(1,:)
  do i = 1, 4
    call dft(Bl(i,:),nn,D,1)
  end do
  deallocate(B1,B2)

  !* normalization
  if (present(Alg))  Alg = 0d0
  if (present(Alc))  Alc = 0d0
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    if (present(Alg)) then
      iAlg = (ll(1,n)**2*(Al(1,n)+Al(4,n)) + ll(1,n)*ll(2,n)*(Al(2,n)+Al(5,n)) + ll(2,n)**2*(Al(3,n)+Al(6,n)))*0.5d0 &
        + ll(1,n)**2*Bl(1,n) + ll(2,n)**2*Bl(2,n) + ll(1,n)*ll(2,n)*sum(Bl(3:4,n))
      if(iAlg>0d0) Alg(n) = 2d0/iAlg
    end if
    if (present(Alc)) then
      iAlc = (ll(2,n)**2*(Al(1,n)+Al(4,n)) - ll(1,n)*ll(2,n)*(Al(2,n)+Al(5,n)) + ll(1,n)**2*(Al(3,n)+Al(6,n)))*0.5d0 &
        + ll(2,n)**2*Bl(1,n) + ll(1,n)**2*Bl(2,n) - ll(1,n)*ll(2,n)*sum(Bl(3:4,n))
      if(iAlc>0d0) Alc(n) = 2d0/iAlc
    end if
  end do

  deallocate(Al)

end subroutine dn0flat_tt


end module recflat


