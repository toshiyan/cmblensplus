!////////////////////////////////////////////////////!
! * Normalization of quadratic lens reconstruction
!////////////////////////////////////////////////////!

module rec_lens
  use alkernel, only: kernels_lens
  use lapack95, only: inv_lapack
  implicit none

  private kernels_lens
  private inv_lapack

contains


subroutine qtt(lmax,rlmin,rlmax,fC,OCT,Ag,Ac)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, G0

  write(*,*) 'norm qTT (lens)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0 / OCT(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  S0 = 0d0
  call kernels_lens(rL,W1,W2,S0,'S0')

  W2 = W1 * fC(rlmin:rlmax)
  G0 = 0d0
  call kernels_lens(rL,W2,W2,G0,'G0')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (S0(1,l)+G0(1,l)/=0d0)  Ag(l) = 1d0/(S0(1,l)+G0(1,l))
    if (S0(2,l)+G0(2,l)/=0d0)  Ac(l) = 1d0/(S0(2,l)+G0(2,l))
  end do
  Ac(1) = 0d0

end subroutine qtt


subroutine qte(lmax,rlmin,rlmax,fC,OCT,OCE,Ag,Ac)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCT, OCE
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, Sp, Gc

  write(*,*) 'norm qTE (lens)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0/OCT(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2/OCE(rlmin:rlmax)
  S0 = 0d0
  call kernels_lens(rL,W1,W2,S0,'S0')

  W1 = 1d0/OCE(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2/OCT(rlmin:rlmax)
  Sp = 0d0
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  W1 = fC(rlmin:rlmax)/OCT(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)/OCE(rlmin:rlmax)
  Gc = 0d0
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (S0(1,l)+Sp(1,l)+2d0*Gc(1,l)/=0d0)  Ag(l) = 1d0/(S0(1,l)+Sp(1,l)+2d0*Gc(1,l))
    if (S0(2,l)+Sp(2,l)+2d0*Gc(2,l)/=0d0)  Ac(l) = 1d0/(S0(2,l)+Sp(2,l)+2d0*Gc(2,l))
  end do
  Ac(1) = 0d0

end subroutine qte


subroutine qtb(lmax,rlmin,rlmax,fC,OCT,OCB,Ag,Ac)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCT, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sm

  write(*,*) 'norm qTB (lens)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2/OCT(rlmin:rlmax)
  Sm = 0d0
  call kernels_lens(rL,W1,W2,Sm,'Sm')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (Sm(1,l)/=0d0)  Ag(l) = 1d0/Sm(1,l)
    if (Sm(2,l)/=0d0)  Ac(l) = 1d0/Sm(2,l)
  end do
  Ag(1) = 0d0

end subroutine qtb


subroutine qee(lmax,rlmin,rlmax,fC,OCE,Ag,Ac)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCE
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sp, Gp

  write(*,*) 'norm qEE (lens)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0/OCE(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  Sp = 0d0
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  W2 = W1 * fC(rlmin:rlmax)
  Gp = 0d0
  call kernels_lens(rL,W2,W2,Gp,'Gp')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (Sp(1,l)+Gp(1,l)/=0d0)  Ag(l) = 1d0/(Sp(1,l)+Gp(1,l))
    if (Sp(2,l)+Gp(2,l)/=0d0)  Ac(l) = 1d0/(Sp(2,l)+Gp(2,l))
  end do
  Ac(1) = 0d0

end subroutine qee


subroutine qeb(lmax,rlmin,rlmax,fC,OCE,OCB,Ag,Ac)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in) , dimension(0:rlmax) :: fC, OCE, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sm

  write(*,*) 'norm qEB (lens)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = fC(rlmin:rlmax)**2 / OCE(rlmin:rlmax)
  Sm = 0d0
  call kernels_lens(rL,W1,W2,Sm,'Sm')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if (Sm(1,l)/=0d0)  Ag(l) = 1d0/Sm(1,l)
    if (Sm(2,l)/=0d0)  Ac(l) = 1d0/Sm(2,l)
  end do
  Ag(1) = 0d0

end subroutine qeb


subroutine qbb(lmax,rlmin,rlmax,fC,OCB,Ag,Ac)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCB
  double precision, intent(out), dimension(0:lmax) :: Ag, Ac
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sp, Gp

  write(*,*) 'norm qBB (lens)'
  rL = (/rlmin,rlmax/)

  W1 = 1d0/OCB(rlmin:rlmax)
  W2 = W1 * fC(rlmin:rlmax)**2
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  W2 = W1 * fC(rlmin:rlmax)
  call kernels_lens(rL,W2,W2,Gp,'Gp')

  Ag = 0d0
  Ac = 0d0
  do l = 1, lmax
    if(Sp(1,l)+Gp(1,l)/=0d0) Ag(l) = 1d0/(Sp(1,l)+Gp(1,l))
    if(Sp(2,l)+Gp(2,l)/=0d0) Ac(l) = 1d0/(Sp(2,l)+Gp(2,l))
  end do
  Ac(1) = 0d0

end subroutine qbb


subroutine qttte(lmax,rlmin,rlmax,fCTT,fCTE,OCT,OCE,OCTE,Ig,Ic)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT,fCTE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: S0, Gc, G0, Sc

  write(*,*) 'norm qTTTE (lens)'

  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  do l = rlmin, rlmax
    W1(l) = 1d0/OCT(l)
    W2(l) = fCTT(l)*fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,S0,'S0')

  do l = rlmin, rlmax
    W1(l) = fCTE(l)/OCT(l)
    W2(l) = fCTT(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)/OCT(l)
  end do
  call kernels_lens(rL,W1,W2,G0,'G0')

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*fCTE(l)/OCT(l)
  end do
  call kernels_lens(rL,W1,W2,Sc,'Sc')

  do l = 1, lmax
    Ig(l) = S0(1,l)+Gc(1,l)+G0(1,l)+Sc(1,l)
    Ic(l) = S0(2,l)+Gc(2,l)+G0(2,l)+Sc(2,l)
  end do

end subroutine qttte


subroutine qttee(lmax,rlmin,rlmax,fCTT,fCEE,OCT,OCE,OCTE,Ig,Ic)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCTT,fCEE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sc, Gc

  write(*,*) 'norm qTTEE (lens)'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0
  W1 = 0d0
  W2 = 0d0

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Sc,'Sc')

  do l = rlmin, rlmax
    W1(l) = fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTT(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  do l = 1, lmax
    Ig(l) = Sc(1,l) + Gc(1,l)
  end do

end subroutine qttee


subroutine qteee(lmax,rlmin,rlmax,fCEE,fCTE,OCT,OCE,OCTE,Ig,Ic)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCTE,OCT,OCE,OCTE
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(2,lmax) :: Sc,Gp,Gc,Sp

  write(*,*) 'norm qTEEE (lens)'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  do l = rlmin, rlmax
    W1(l) = OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTE(l)*fCEE(l)/OCE(l)
  end do
  call kernels_lens(rL,W1,W2,Sc,'Sc')

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCEE(l)/OCE(l)
  end do
  call kernels_lens(rL,W1,W2,Gp,'Gp')

  do l = rlmin, rlmax
    W1(l) = fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCTE(l)/OCE(l)
  end do
  call kernels_lens(rL,W1,W2,Gc,'Gc')

  do l = rlmin, rlmax
    W1(l) = 1d0/OCE(l)
    W2(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Sp,'Sp')

  do l = 1, lmax
    Ig(l) = (Sc(1,l)+Gp(1,l)+Gc(1,l)+Sp(1,l))
    Ic(l) = (Sc(2,l)+Gp(2,l)+Gc(2,l)+Sp(2,l))
  end do

end subroutine qteee


subroutine qtbeb(lmax,rlmin,rlmax,fCEE,fCBB,fCTE,OCT,OCE,OCB,OCTE,Ig,Ic)
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE,fCBB,fCTE,OCT,OCE,OCTE,OCB
  double precision, intent(out), dimension(0:lmax) :: Ig, Ic
  !internal
  integer :: l, rL(2)
  double precision, dimension(2,lmax) :: Gm, Sm
  double precision, dimension(rlmin:rlmax) :: W1, W2

  write(*,*) 'norm qTBEB (lens)'
  rL = (/rlmin,rlmax/)
  Ig = 0d0
  Ic = 0d0

  do l = rlmin, rlmax
    W1(l) = fCTE(l)*OCTE(l)/(OCT(l)*OCE(l))
    W2(l) = fCBB(l)/OCB(l)
  end do
  call kernels_lens(rL,W1,W2,Gm,'Gm')

  do l = rlmin, rlmax
    W1(l) = 1d0/OCB(l)
    W2(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCT(l)*OCE(l))
  end do
  call kernels_lens(rL,W1,W2,Sm,'Sm')

  Ig(1:lmax) = Sm(1,1:lmax) + Gm(1,1:lmax)
  Ic(1:lmax) = Sm(2,1:lmax) + Gm(2,1:lmax)

end subroutine qtbeb


subroutine qmv(lmax,QDO,Al,Il,MV)
!* Compute MV estimator from pre-computed normalization (Al) and correlations (Il)
! Oder of QDO: TT, TE, TB, EE, EB, BB
  implicit none
  ![input]
  ! lmax         --- Maximum multipole of MV(l)
  ! QDO(6)       --- T/F array, specifying which estimators to be used for MV
  ! Al(6,lmax+1) --- Normalization
  ! Il(4,lmax+1) --- Correlations
  integer, intent(in) :: lmax
  logical, intent(in), dimension(6) :: QDO
  double precision, intent(in), dimension(6,0:lmax) :: Al
  double precision, intent(in), dimension(4,0:lmax) :: Il
  double precision, intent(out), dimension(0:lmax) :: MV
  !internal
  integer :: QTT = 1, QTE = 2, QTB = 3, QEE = 4, QEB = 5, QBB = 6
  integer :: X, Y, qmax, i, id(6), l
  double precision, allocatable :: M(:,:)

  write(*,*) 'norm qMV (lens)'
  id = 0

  !* set ids
  do X = 1, 6
    if(QDO(X)) id(X) = 1 + maxval(id)
  end do
  qmax = maxval(id)

  MV = 0d0

  do l = 1, lmax

    !* noise covariance
    allocate(M(qmax,qmax));  M = 0d0

    if(QDO(QTT).and.QDO(QTE)) M(id(QTT),id(QTE)) = Il(1,l)*Al(QTT,l)*Al(QTE,l)
    if(QDO(QTT).and.QDO(QEE)) M(id(QTT),id(QEE)) = Il(2,l)*Al(QTT,l)*Al(QEE,l)
    if(QDO(QTE).and.QDO(QEE)) M(id(QTE),id(QEE)) = Il(3,l)*Al(QTE,l)*Al(QEE,l)
    if(QDO(QTB).and.QDO(QEB)) M(id(QTB),id(QEB)) = Il(4,l)*Al(QTB,l)*Al(QEB,l)

    do X = 1, 6
      if (QDO(X)) M(id(X),id(X)) = Al(X,l)
      do Y = X + 1, 6
        if(QDO(X).and.QDO(Y)) M(id(Y),id(X)) = M(id(X),id(Y))
      end do
    end do 
    call inv_lapack(M)

    MV(l) = 1d0/sum(M)

    deallocate(M)

  end do

end subroutine qmv


end module rec_lens

