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


