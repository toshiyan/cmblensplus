!////////////////////////////////////////////////////!
! * All and iterative estimator
!////////////////////////////////////////////////////!

module nldd_interface
  use myconst, only: pi, TT, TE, EE, BB, dd
  !use mylapack95, only: inv_lapack
  use nldd_lens, only: AlTT, AlTE, AlTB, AlEE, AlEB, AlBB, ILTTTE, ILTTEE, ILTEEE, ILTBEB
  use nldd_delens, only: AlEB_iter
  implicit none

  private pi, TT, TE, EE, BB, dd
  private inv_lapack
  private AlTT, AlTE, AlTB, AlEE, AlEB, AlBB, AlEB_iter, ILTTTE, ILTTEE, ILTEEE, ILTBEB

  !* public variables
  integer, parameter :: QTT = 1, QTE = 2, QTB = 3, QEE = 4, QEB = 5, QBB = 6, QMV = 7
  integer, parameter :: QTTTE = 1, QTTEE = 2, QTEEE = 3, QTBEB = 4

contains

#ifdef all
subroutine OPTCOMB(QDO,Al,Il,Nl)
  implicit none
  !I/O
  logical, intent(in) :: QDO(:)
  double precision, intent(in) :: Il(:)
  double precision, intent(inout) :: Al(:)
  double precision, intent(out), optional :: Nl(:)
  !internal
  integer :: X, Y, qmax, i, id(6)
  double precision, dimension(:,:), allocatable :: M

  id = 0

  !* set ids
  do X = 1, 6
    if(QDO(X)) id(X) = 1 + maxval(id)
  end do
  qmax = maxval(id)

  !* noise covariance
  allocate(M(qmax,qmax));  M = 0d0

  if(QDO(QTT).and.QDO(QTE)) M(id(QTT),id(QTE)) = Il(QTTTE)*Al(QTT)*Al(QTE)
  if(QDO(QTT).and.QDO(QEE)) M(id(QTT),id(QEE)) = Il(QTTEE)*Al(QTT)*Al(QEE)
  if(QDO(QTE).and.QDO(QEE)) M(id(QTE),id(QEE)) = Il(QTEEE)*Al(QTE)*Al(QEE)
  if(QDO(QTB).and.QDO(QEB)) M(id(QTB),id(QEB)) = Il(QTBEB)*Al(QTB)*Al(QEB)

  do X = 1, 6
    if (QDO(X)) M(id(X),id(X)) = Al(X)
    do Y = X + 1, 6
      if(QDO(X).and.QDO(Y)) M(id(Y),id(X)) = M(id(X),id(Y))
    end do
  end do 
  call INV_LAPACK(M)

  Al(QMV) = 1d0/sum(M)
  if (present(Nl)) Nl = sum(M,dim=2)
  deallocate(M)

end subroutine OPTCOMB
#endif all

subroutine AL_INTERFACE(rL,dL,fC,OC,Alg,Alc,QDO,itern)
  implicit none
  !I/O
  logical, intent(in) :: QDO(1:7)
  integer, intent(in) :: rL(1:2), dL(1:2)
  integer, intent(in), optional :: itern
  double precision, intent(out) :: Alg(:,:), Alc(:,:)
  double precision, intent(in), dimension(:,:) :: fC, OC
  !internal
  integer :: i, n, l
  double precision :: conv, ratio
  double precision, dimension(:), allocatable :: AlgEB, rCBB
  double precision, dimension(:,:), allocatable :: Ilg, Ilc

  conv = 1d-6

  !//// interface ////!
  ! TT
  if (QDO(QTT))  call ALTT(rL,dL,Alg(QTT,:),Alc(QTT,:),fC(TT,:),OC(TT,:))
  ! TE
  if (QDO(QTE))  call ALTE(rL,dL,Alg(QTE,:),Alc(QTE,:),fC(TE,:),OC(TT,:),OC(EE,:))
  ! EE
  if (QDO(QEE))  call ALEE(rL,dL,Alg(QEE,:),Alc(QEE,:),fC(EE,:),OC(EE,:))
  ! TB
  if (QDO(QTB))  call ALTB(rL,dL,Alg(QTB,:),Alc(QTB,:),fC(TE,:),OC(TT,:),OC(BB,:))
  ! BB
  if (QDO(QBB))  call ALBB(rL,dL,Alg(QBB,:),Alc(QBB,:),fC(BB,:),OC(BB,:))
  ! EB
  if (QDO(QEB)) then
    if (present(itern).and.itern>0) then
      call AlEB_iter(itern,rL,dL,fC(EE,:),fC(dd,:),OC(EE,:),OC(BB,:),Alg(QEB,:),Alc(QEB,:))
    else
      call ALEB(rL,dL,Alg(QEB,:),Alc(QEB,:),fC(EE,:),OC(EE,:),OC(BB,:))
    end if
  end if

  !MV
  if(QDO(QMV)) then
    allocate(Ilg(4,dL(2)),Ilc(4,dL(2)))
    ! TT x TE
    if(QDO(QTT).and.QDO(QTE)) call ILTTTE(rL,dL,Ilg(QTTTE,:),Ilc(QTTTE,:),fC(TT,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:))
    ! TT x EE
    if(QDO(QTT).and.QDO(QEE)) call ILTTEE(rL,dL,Ilg(QTTEE,:),Ilc(QTTEE,:),fC(TT,:),fC(EE,:),OC(TT,:),OC(EE,:),OC(TE,:))
    ! TE x EE
    if(QDO(QTE).and.QDO(QEE)) call ILTEEE(rL,dL,Ilg(QTEEE,:),Ilc(QTEEE,:),fC(EE,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:))
    ! TB x EB
    if(QDO(QTB).and.QDO(QEB)) call ILTBEB(rL,dL,Ilg(QTBEB,:),Ilc(QTBEB,:),fC(EE,:),fC(BB,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(BB,:),OC(TE,:))
    do l = dL(1), dL(2)
      call OPTCOMB(QDO,Alg(:,l),Ilg(:,l))
      call OPTCOMB(QDO,Alc(:,l),Ilc(:,l))
    end do
    deallocate(Ilg,Ilc)
  end if

end subroutine AL_INTERFACE

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
  double complex, allocatable :: Al(:,:), Bl(:,:), A1(:,:), A2(:,:), B1(:,:), B2(:,:)

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






end module nldd_interface

