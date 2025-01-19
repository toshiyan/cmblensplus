!////////////////////////////////////////////////////!
! * Delensing 
!////////////////////////////////////////////////////!

module delensing
  use constants, only: pi
  use funcs,   only: WIGd_INI!, WIGD_RECURSION
  use general, only: gauss_legendre_params, gl_initialize, gl_finalize, savetxt, linspace
  use alkernel, only: conv_egrad, zeta, zeta_l
  implicit none

  !local parameters
  double precision, parameter :: cx(1:2) = [1d0,-1d0], fourpi = 4d0*pi
  private cx, fourpi

  private pi
  private wigd_ini
  private gauss_legendre_params, gl_initialize, gl_finalize, savetxt, linspace
  private conv_egrad, zeta, zeta_l

contains


!subroutine RES_CLBB(oL,dLE,dLP,CE,Cp,CB,WE,Wp,NE,Np)
!* residual ClBB = ClBB^lin - ClBB^est
!
!  implicit none
! [inputs]  
!   oL  --- multipole range of residual ClBB
!   dLE --- multipole range of E modes
!   dLP --- multipole range of lensing potential
!   CE, Cp --- power spectrum of E-mode and lensing pontential
!  integer, intent(in) :: oL(2), dLE(2), dLP(2)
!  double precision, intent(in), dimension(:) :: CE, Cp
!
! (optional)
!   WE, Wp --- Wiener filters of E-mode and lensing potential
!   NE, Np --- Noise power spectra of E-mode and lensing potential
!  double precision, intent(in), dimension(:), optional :: WE, Wp, NE, Np
!
! [outputs]
!   CB --- residual B-mode
!  double precision, intent(out) :: CB(:)
!
! [internal]
!  integer :: l, lmin, lmax
!  double precision :: CBW(oL(2)), CBL(oL(2))
!  double precision, dimension(dLE(1):dLE(2)) :: A
!  double precision, dimension(dLP(1):dLP(2)) :: B

!  CBW = 0d0

!  if (present(WE).and.present(Wp)) then
!    A = CE(dLE(1):dLE(2)) * WE(dLE(1):dLE(2))
!    B = Cp(dLP(1):dLP(2)) * Wp(dLP(1):dLP(2))
!    call conv_egrad(oL,dLE,dLP,A,B,CBW)
!  end if

!  if (present(NE).and.present(Np)) then
!    do l = dLE(1), dLE(2)
!      A(l) = CE(l)**2/(CE(l)+NE(l))
!    end do
!    do l = dLP(1), dLP(2)
!      B(l) = Cp(l)**2/(Cp(l)+Np(l))
!    end do
!    call conv_egrad(oL,dLE,dLP,A,B,CBW)
!  end if

! Lensing B-mode power spectrum
!  call conv_egrad(oL,(/1,dLE(2)/),(/1,dLP(2)/),CE(1:dLE(2)),Cp(1:dLP(2)),CBL)

!  CB = CBL
!  CB(oL(1):oL(2)) = CBL(oL(1):oL(2)) - CBW(oL(1):oL(2))

!end subroutine RES_CLBB


!subroutine CLBB_LIN(eL,dLE,dLP,CE,CP,CB)
! * Lensing B-mode power spectrum as a convolution of ClEE and Clpp
!  implicit none
!
! [input]
!   eL --- multipole range of residual ClBB
!   dLE, dLP --- multipole range of E modes and lensing potential
!   CE, CP --- power spectrum of E-mode and lensing pontential
!  integer, intent(in) :: eL(2), dLE(2), dLP(2)
!  double precision, intent(in), dimension(:) :: CE, CP
!
! [output]
!   CB --- residual B-mode
!  double precision, intent(out) :: CB(:)

!  call conv_egrad(eL,dLE,dLP,CE(dLE(1):dLE(2)),CP(dLP(1):dLP(2)),CB)

!end subroutine CLBB_LIN


subroutine CLBB_EST(eL,dLE,dLP,CE,Cp,NE,Np,Cl)
! * Estimate of lensing B-mode power spectrum (Noise power spectra as inputs)
  implicit none
!
! [input]
!   eL --- multipole range of residual ClBB
!   dL --- multipole range of delensing
!   CE --- E-mode power spectrum
!   Cp --- lensing potential power spectrum
!   NE --- E-mode noise power spectrum
!   Np --- lensing potential noise power spectrum
  integer, intent(in) :: dLE(2), dLP(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, Cp, NE, Np
!
! [output]
!   Cl --- estimated lensing B-mode power spectrum
  double precision, intent(out) :: Cl(:)
!
! [internal]
  integer :: l
  double precision, dimension(dLE(1):dLE(2)) :: WE
  double precision, dimension(dLP(1):dLP(2)) :: WP

  do l = dLE(1), dLE(2)
    WE(l) = CE(l)**2/(CE(l)+NE(l))
  end do
  do l = dLP(1), dLP(2)
    Wp(l) = Cp(l)**2/(Cp(l)+Np(l))
  end do

  call conv_egrad(eL,dLE,dLP,WE,Wp,Cl)

end subroutine CLBB_EST


subroutine Delensing_Bias(eL,dL,LE,LB,Cp,DBL,NP1,Ag,NP2)
!* Compute the dominant terms of the delensing bias
  implicit none
  integer, intent(in) :: dL(2), eL(2)
  double precision, intent(in), dimension(:) :: LE,LB,Cp,NP1,Ag
  double precision, intent(in), optional :: NP2(:)
  double precision, intent(out) :: DBL(:)
  !internal
  integer :: l
  double precision, dimension(eL(2)) :: wBB, DLs

  !* delensing bias (Dl and sum of bias terms)

  DLs = 0d0
  DBl = 0d0

  call DLBB(eL,dL,LE,LE,LE+NP1,LB+NP1,Cp,Ag,DLS,LE+NP2)
  call CLBB_EST(eL,dL,dL,LE,Cp,NP2,Ag,wBB)

  do l = eL(1), eL(2)
    DBl(l) = DLS(l)*(-2d0*LB(l)+2d0*wBB(l)+DLS(l)*(LB(l)+NP1(l)))
  end do

end subroutine Delensing_Bias


subroutine DlBB(eL,dL,CE,LE,OE1,OB1,Cp,Np,Cl,OE2)
!* Computing a power spectrum (Dl) where B^bias_lm ~ D_l B_lm
  implicit none
!
! [input]
!   eL --- multipole range of residual ClBB
!   dL --- multipole range of delensing
!   CE --- E-mode power spectrum
!   Cp --- lensing potential power spectrum
!   NE --- E-mode noise power spectrum
!   Np --- lensing potential noise power spectrum
  integer, intent(in) :: dL(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, LE, OE1, OB1, Cp, Np
!
! (optional)
  double precision, intent(in), dimension(:), optional :: OE2
!
! [output]
!   Cl --- delensing bias
  double precision, intent(out) :: Cl(:)
!
! [internal]
  integer :: l
  double precision, dimension(dL(1):dL(2)) :: WE, Wp

  do l = dL(1), dL(2)
    WE(l) = 0d0
    if(l>=dL(1)) WE(l) = (LE(l)/OE1(l))*CE(l)
    if(present(OE2)) WE(l) = WE(l)*LE(l)/OE2(l)
    Wp(l) = Cp(l)*Np(l)/(Cp(l)+Np(l))
  end do
  call conv_egrad(eL,dL,dL,WE,Wp,Cl)

  Cl(eL(1):eL(2)) = Cl(eL(1):eL(2))/OB1(eL(1):eL(2))

end subroutine DlBB


!//// covariance of lensed/delensed ClBB ////!

subroutine dCBdCE(oL,eL,Cp,X)
! * Derivative of B-mode Cl in terms of the E-mode Cl
!
  implicit none
! [input]  
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   Cp --- power spectrum of the lensing pontential
  integer, intent(in) :: oL(1:2), eL(1:2)
  double precision, intent(in), dimension(eL(1):eL(2)) :: Cp
!
! [output]
!   X --- derivative of ClBB in terms of Clpp
  double precision, intent(out) :: X(eL(2),oL(2))
!
! [internal]
  type(gauss_legendre_params) :: GL
  integer :: i, l
  double precision :: mu, al, c1_inv, c2p, c2m, c3
  double precision, dimension(2) :: ZP, d22_sup, d22_mid, d22_inf
  double precision, dimension(eL(2)) :: Ip, Im
  double precision, dimension(eL(2),2) :: ZE33, ZE31, ZE11
  double precision, dimension(eL(1):eL(2)) :: al0, alm, alp

  do l = eL(1), eL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp(l) = dsqrt((al-2)*(al+3))
    alm(l) = dsqrt((al+2)*(al-1))
  end do

  !* Gauss-Legendre Integration
  call gl_initialize(GL,int((3*eL(2)+1)/2),1d-15)
  X = 0d0
  do i = 1, GL%n
    mu = GL%z(i)
    call ZETA_L(3,3,eL,alp**2,mu,ZE33)
    call ZETA_L(3,1,eL,alp*alm,mu,ZE31)
    call ZETA_L(1,1,eL,alm**2,mu,ZE11)
    call ZETA(1,1,eL,Cp*al0**2,mu,ZP)
    Ip = (ZE33(:,1)*ZP(1)+2*ZE31(:,1)*ZP(2)+ZE11(:,1)*ZP(1))
    Im = (ZE33(:,2)*ZP(2)+2*ZE31(:,2)*ZP(1)+ZE11(:,2)*ZP(2))
    d22_sup = 0d0 ;  d22_mid = 0d0 ;  d22_inf = 0d0
    do l = 2, oL(2) !recursion of wigner d function
      if(l==2) then 
        d22_sup = wigd_ini(2,2,mu)
      else
        al = dble(l)
        c1_inv = (al*(2d0*al-1d0))/(al**2-4d0)
        c2p = 4d0/((al-1d0)*al) - mu
        c2m = -4d0/((al-1d0)*al) - mu
        c3 = ((al-1d0)**2-4d0)/((al-1d0)*(2d0*al-1d0))
        d22_sup(1) = -(c2p*d22_mid(1) + c3*d22_inf(1))*c1_inv
        d22_sup(2) = -(c2m*d22_mid(2) + c3*d22_inf(2))*c1_inv
      end if
      X(:,l) = X(:,l) + (Ip(:)*d22_sup(1)-Im(:)*d22_sup(2))*GL%w(i)*pi*0.25d0
      d22_inf = d22_mid
      d22_mid = d22_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine dCBdCE


subroutine dCBdCp(oL,eL,CE,X)
! * derivative of B-mode Cl in terms of the lensing potential Cl
!
  implicit none
! [inputs]  
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   CE --- power spectrum of the E-modes
  integer, intent(in) :: oL(1:2), eL(1:2)
  double precision, intent(in), dimension(eL(1):eL(2)) :: CE
!
! [outputs]
!   X --- derivative of ClBB in terms of Clpp
  double precision, intent(out) :: X(eL(2),oL(2))
!
! [internal]
  type(gauss_legendre_params) :: GL
  integer :: i, l
  double precision :: mu, al, c1_inv, c2p, c2m, c3
  double precision, dimension(2) :: d22_sup, d22_mid, d22_inf, ZE33, ZE31, ZE11
  double precision, dimension(eL(2)) :: Ip, Im
  double precision, dimension(eL(2),2) :: ZP
  double precision, dimension(eL(1):eL(2)) :: al0, alm, alp

  do l = eL(1), eL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp(l) = dsqrt((al-2)*(al+3))
    alm(l) = dsqrt((al+2)*(al-1))
  end do

  X = 0d0
  !* Gauss-Legendre Integration
  call gl_initialize(GL,int((3*eL(2)+1)/2),1d-15)

  do i = 1, GL%n
    mu = GL%z(i)
    call ZETA(3,3,eL,CE*alp**2,mu,ZE33)
    call ZETA(3,1,eL,CE*alp*alm,mu,ZE31)
    call ZETA(1,1,eL,CE*alm**2,mu,ZE11)
    call ZETA_L(1,1,eL,al0**2,mu,ZP)
    Ip = (ZE33(1)*ZP(:,1)+2*ZE31(1)*ZP(:,2)+ZE11(1)*ZP(:,1))
    Im = (ZE33(2)*ZP(:,2)+2*ZE31(2)*ZP(:,1)+ZE11(2)*ZP(:,2))
    d22_sup = 0d0 ;  d22_mid = 0d0 ;  d22_inf = 0d0
    do l = 2, oL(2) !recursion of wigner d function
      if(l==2) then 
        d22_sup = wigd_ini(2,2,mu)
      else
        al = dble(l)
        c1_inv = (al*(2d0*al-1d0))/(al**2-4d0)
        c2p = 4d0/((al-1d0)*al) - mu
        c2m = -4d0/((al-1d0)*al) - mu
        c3 = ((al-1d0)**2-4d0)/((al-1d0)*(2d0*al-1d0))
        d22_sup(1) = -(c2p*d22_mid(1) + c3*d22_inf(1))*c1_inv
        d22_sup(2) = -(c2m*d22_mid(2) + c3*d22_inf(2))*c1_inv
      end if
      X(:,l) = X(:,l) + (Ip(:)*d22_sup(1)-Im(:)*d22_sup(2))*GL%w(i)*pi*0.25d0
      d22_inf = d22_mid
      d22_mid = d22_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine dCBdCp


end module delensing


