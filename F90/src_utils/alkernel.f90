!////////////////////////////////////////////////////!
! * Kernel Functions of Normalization
!////////////////////////////////////////////////////!

module alkernel
  use constants, only: pi
  use funcs,   only: WIGd_INI!, WIGD_RECURSION
  use general, only: gauss_legendre_params, gl_initialize, gl_finalize, savetxt, linspace
  implicit none

  !local parameters
  double precision, parameter :: cx(1:2) = [1d0,-1d0], fourpi = 4d0*pi
  private cx, fourpi

  private pi
  private wigd_ini
  private gauss_legendre_params, gl_initialize, gl_finalize, savetxt, linspace

contains


subroutine get_lfac(lmax,lfac,lk2)
  implicit none
  !I/O
  double precision, intent(out), dimension(lmax) :: lk2
  integer, intent(in) :: lmax
  character(1), intent(in) :: lfac
  integer :: l

  lk2 = 1d0
  if (lfac=='k') then
    do l = 1, lmax
      lk2(l) = (dble(l*(l+1))/2d0)**2
    end do
  end if

end subroutine get_lfac


subroutine Zeta(m1,m2,rL,A,mu,Z)
! * Computing Zeta_{m1,m2}(A,\mu) = \sum_{l=m1}^{lmax} A(l) * (2l+1)*d_{m1m2}^l/4pi
! - return Zeta(|m1|,|m2|) and Zeta(|m1|,-|m2|)
! assuming |m1|>=|m2|
  implicit none
  !I/O
  integer, intent(in) :: m1, m2, rL(1:2)
  double precision, intent(in) :: mu, A(rL(1):rL(2))
  double precision, intent(out) :: Z(1:2)
  !internal
  integer :: l
  double precision, dimension(1:2) :: dmm, d_mid, d_inf
  double precision :: c1_inv,c3,c2p,c2m,x
  double precision :: aL,aL2,aLL2,aM11,aM22,aM12

  !* initialize
  dmm = 0d0 ;  d_mid = 0d0 ;  d_inf = 0d0 
  aM11 = dble(m1)**2
  aM22 = dble(m2)**2
  aM12 = dble(m1)*dble(m2)

  !* l = m1
  dmm = wigd_ini(m1,m2,mu)
  Z = 0d0
  if (m1>=rL(1)) Z = Z + A(m1)*dble(2*m1+1)/fourpi*dmm
  d_inf = d_mid
  d_mid = dmm

  !* l > m1
  do l = m1+1, rL(2) ! recusrion formula of the wigner d matrix
    aL = dble(l)
    aL2 = aL**2
    aLL2 = (aL-1d0)**2
    c1_inv = (2d0*aL2-aL)/dsqrt((aL2-aM11)*(aL2-aM22))
    if (l==m1+1) then 
      c3 = 0d0
      if (m1==0) then
        c2p = - mu
        c2m = - mu
      else
        x = aM12/(aL2-aL)
        c2p = x - mu
        c2m = -x - mu
      end if
    else
      x = aM12/(aL2-aL)
      c2p = x - mu
      c2m = - x - mu
      c3 = dsqrt((aLL2-aM11)*(aLL2-aM22))/((aL-1d0)*(2d0*aL-1d0))
    end if
    dmm(1) = -(c2p*d_mid(1) + c3*d_inf(1))*c1_inv
    dmm(2) = -(c2m*d_mid(2) + c3*d_inf(2))*c1_inv
    if(l>=rL(1))  Z = Z + A(l)*(2d0*aL+1d0)/fourpi*dmm
    d_inf = d_mid
    d_mid = dmm
  end do

end subroutine Zeta


subroutine ZETA_L(m1,m2,rL,A,mu,Z)
! * Zeta at each multipoles
!     Zeta_{m1,m2}(A,\mu,L) = A(l) * (2l+1)*d_{m1m2}^l/4pi
!
  implicit none
  !I/O
  integer, intent(in) :: m1, m2, rL(1:2)
  double precision, intent(in) :: mu, A(rL(1):rL(2))
  double precision, intent(out) :: Z(rL(2),1:2)
  !internal
  integer :: l
  double precision, dimension(1:2) :: dmm, d_mid, d_inf
  double precision :: c1_inv,c3,c2p,c2m,x
  double precision :: aL,aL2,aLL2,aM11,aM22,aM12

  !* initialize
  dmm = 0d0 ;  d_mid = 0d0 ;  d_inf = 0d0 
  aM11 = dble(M1)**2
  aM22 = dble(M2)**2
  aM12 = dble(M1)*dble(M2)

  !* l = m1
  dmm = wigd_ini(m1,m2,mu)
  Z = 0d0
  if (m1>=rL(1)) Z(m1,:) = A(m1)*dble(2*m1+1)/fourpi*dmm
  d_inf = 0d0
  d_mid = dmm

  !* l > m1
  do l = m1+1, rL(2) ! recusrion formula of the wigner d matrix
    aL = dble(L)
    aL2 = aL**2
    aLL2 = (aL-1d0)**2
    c1_inv = (2d0*aL2-aL)/dsqrt((aL2-aM11)*(aL2-aM22))
    if(l==m1+1) then 
      c3 = 0d0
      if (m1==0) then
        c2p = - mu
        c2m = - mu
      else
        x = aM12/(aL2-aL)
        c2p = x - mu
        c2m = - x - mu
      end if
    else
      x = aM12/(aL2-aL)
      c2p = x - mu
      c2m = - x - mu
      c3 = dsqrt((aLL2-aM11)*(aLL2-aM22))/((aL-1d0)*(2d0*aL-1d0))
    end if
    dmm(1) = -(c2p*d_mid(1) + c3*d_inf(1))*c1_inv
    dmm(2) = -(c2m*d_mid(2) + c3*d_inf(2))*c1_inv
    if(l>=rL(1))  Z(l,:) = A(l)*(2d0*aL+1d0)/fourpi*dmm
    d_inf = d_mid
    d_mid = dmm
  end do

end subroutine ZETA_L


!/////////////////////////////////////////////////////////////////////////////!
! Delensing

subroutine conv_egrad(oL,dLE,dLP,WE,WP,X)
! * Kernel of lensing B-mode power spectrum
!
! The Kernel function is given by
!
!   (1/(2l+1)) \sum_{L1L2} (S^-_{lL1L2})^2 WE_L1 WP_L2
!
! see Eqs.(2.7) and (A.11) in Namikawa & Nagata 2014
!
! inputs  : oL --- multipole range of X
!         : dLE, dLP --- multipole range of WE and WP
!         : WE, WP --- functions to be convolved
! outputs : X --- kernel function
  implicit none
  !I/O
  integer, intent(in) :: oL(1:2), dLE(1:2), dLP(1:2)
  double precision, intent(in), dimension(dLE(1):dLE(2)) :: WE
  double precision, intent(in), dimension(dLP(1):dLP(2)) :: WP
  double precision, intent(out) :: X(:)
  !internal
  type(gauss_legendre_params) :: GL
  integer :: i, l, lmin, lmax
  double precision :: mu, al, Ip, Im, c1_inv, c2p, c2m, c3
  double precision, dimension(2) :: ZE33, ZE31, ZE11, ZP, d22_sup, d22_mid, d22_inf
  double precision, dimension(min(dLE(1),dLP(1)):max(dLE(2),dLP(2))) :: al0, alm, alp

  lmin = min(dLE(1),dLP(1))
  lmax = max(dLE(2),dLP(2))

  do l = lmin, lmax
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp(l) = dsqrt((al-2)*(al+3))
    alm(l) = dsqrt((al+2)*(al-1))
  end do

  X = 0d0

  !* GL quadrature
  call gl_initialize(GL,int((3*lmax+1)/2),1d-15)

  do i = 1, GL%n
    mu = GL%z(i)
    call ZETA(3,3,dLE,WE*alp(dLE(1):dLE(2))**2,mu,ZE33)
    call ZETA(3,1,dLE,WE*alp(dLE(1):dLE(2))*alm(dLE(1):dLE(2)),mu,ZE31)
    call ZETA(1,1,dLE,WE*alm(dLE(1):dLE(2))**2,mu,ZE11)
    call ZETA(1,1,dLP,WP*al0(dLP(1):dLP(2))**2,mu,ZP)
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
      Ip = (ZE33(1)*ZP(1)+2*ZE31(1)*ZP(2)+ZE11(1)*ZP(1))
      Im = (ZE33(2)*ZP(2)+2*ZE31(2)*ZP(1)+ZE11(2)*ZP(2))
      X(l) = X(l) + (Ip*d22_sup(1)-Im*d22_sup(2))*GL%w(i)*pi*0.25d0
      d22_inf = d22_mid
      d22_mid = d22_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine conv_egrad


subroutine Kernels_Lens(rL,WA,WB,X,kernel)
!Kernels of lensing reconstruction noise
! for Gamma, this computes px Gamma = cx**2 Gamma.
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(2)
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:,:)
  !internal
  type(gauss_legendre_params) :: GL
  character(LEN=4) :: selec
  integer :: i, l, lmax
  double precision :: pm, mu, Ip, Im, al, c1_inv, c2p, c2m, c3
  double precision, dimension(rL(1):rL(2)) :: al0, alp2, alm2
  double precision, dimension(2) :: d11, d11_mid, d11_inf, C
  double precision, dimension(2) :: d00, d00_mid, d00_inf
  double precision, dimension(2) :: d10, d10_mid, d10_inf
  double precision, dimension(2) :: ZA00, ZA10, ZA20, ZA21, ZA22, ZA30, ZA32
  double precision, dimension(2) :: ZB00, ZB10, ZB11, ZB31, ZB32, ZB33, ZB21

  !* initialize
  lmax = size(X,dim=2)

  if (kernel=='Sp'.or.kernel=='Gp') pm = 1d0
  if (kernel=='Sm'.or.kernel=='Gm') pm = -1d0

  selec = kernel
  if (kernel=='Sp'.or.kernel=='Sm') selec = 'Spm'
  if (kernel=='Gp'.or.kernel=='Gm') selec = 'Gpm'

  X = 0d0

  do l = rL(1), rL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp2(l) = dsqrt((al-2)*(al+3))
    alm2(l) = dsqrt((al+2)*(al-1))
  end do

  !* GL quadrature
  call gl_initialize(GL,int((3*max(rL(2),lmax)+1)/2),1d-15)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(selec)
    case ('S0')
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(1,1,rL,WB*al0**2,mu,ZB11)
      Ip = ZA00(1)*ZB11(1)
      Im = ZA00(1)*ZB11(2)
    case ('G0')
      call ZETA(1,0,rL,WA*al0,mu,ZA10)
      call ZETA(1,0,rL,WB*al0,mu,ZB10)
      Ip = -ZA10(1)*ZB10(1)
      Im = -Ip
    case ('Sc')
      call Zeta(2,0,rL,WA,mu,ZA20)
      call Zeta(1,1,rL,WB*al0*alm2,mu,ZB11)
      call Zeta(3,1,rL,WB*al0*alp2,mu,ZB31)
      Ip = ZA20(1)*(ZB11(2)+ZB31(1))*0.5d0
      Im = ZA20(1)*(ZB31(2)+ZB11(1))*0.5d0
    case ('Gc')
      call Zeta(3,0,rL,WA*alp2,mu,ZA30)
      call Zeta(1,0,rL,WA*alm2,mu,ZA10)
      call Zeta(2,1,rL,WB*al0,mu,ZB21)
      Ip = -(ZA30(1)*ZB21(2)+ZA10(1)*ZB21(1))*0.5d0
      Im = -(ZA30(1)*ZB21(1)+ZA10(1)*ZB21(2))*0.5d0
    case ('Spm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(3,3,rL,WB*alp2**2,mu,ZB33)
      call Zeta(3,1,rL,WB*alp2*alm2,mu,ZB31)
      call Zeta(1,1,rL,WB*alm2**2,mu,ZB11)
      Ip = ((ZB33(1)+ZB11(1))*ZA22(1)+pm*2d0*ZB31(2)*ZA22(2))/4d0
      Im = (pm*(ZB33(2)+ZB11(2))*ZA22(2)+2d0*ZB31(1)*ZA22(1))/4d0
    case ('Gpm')
      call Zeta(3,2,rL,WA*alp2,mu,ZA32)
      call Zeta(2,1,rL,WA*alm2,mu,ZA21)
      call Zeta(3,2,rL,WB*alp2,mu,ZB32)
      call Zeta(2,1,rL,WB*alm2,mu,ZB21)
      Ip = -(ZA21(1)*ZB32(1)+ZA32(1)*ZB21(1)+pm*(ZA32(2)*ZB32(2)+ZA21(2)*ZB21(2)))/4d0
      Im = (ZA32(1)*ZB32(1)+ZA21(1)*ZB21(1)-pm*(ZA21(2)*ZB32(2)+ZA32(2)*ZB21(2)))/4d0
    case default
      stop 'error: no kernel'
    end select
    do l = 1, lmax ! loop for l
      al = dble(l)
      if (l==1) then 
        !d00_mid  = 1d0 !d^0_00(mu)
        !d10_mid  = 0d0 !d^0_10(mu)
        d11_mid  = 0d0 !d^0_11(mu)
        !d00(1:2) = mu  !d^1_00(mu)
        !d10 = -dsqrt(1-mu**2)/dsqrt(2.d0) !d^1_10(mu)
        d11(1) = (1+mu)*0.5d0 !d^1_11(mu)
        d11(2) = (1-mu)*0.5d0 !d^1_1,-1(mu)
      else
        !* M1=M2=0
        !c1_inv = (2d0*al-1d0)/al
        !c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        !d00(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        !d00(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
        !* M1=M2=1
        c1_inv = al*(2d0*al-1d0)/(al**2-1d0)
        c2p = 1d0/((aL-1d0)*aL) - mu
        c2m = -1d0/((aL-1d0)*aL) - mu
        c3 = ((al-1d0)**2-1d0)/((al-1d0)*(2d0*al-1d0))
        d11(1) = -(c2p*d11_mid(1) + c3*d11_inf(1))*c1_inv
        d11(2) = -(c2m*d11_mid(2) + c3*d11_inf(2))*c1_inv
        !* M1=1,M2=0
        !c1_inv = al*(2d0*al-1d0)/(dsqrt((al**2-1d0))*al)
        !c3 = dsqrt(((al-1d0)**2-1d0))*(al-1d0)/((al-1d0)*(2d0*al-1d0))
        !d10(1) = (mu*d10_mid(1) - c3*d10_inf(1))*c1_inv
        !d10(2) = (mu*d10_mid(2) - c3*d10_inf(2))*c1_inv
      end if
      C(:) = (Ip*d11(1)+cx(:)*Im*d11(2))*dble(l*(l+1))
      X(1:2,l) = X(1:2,l) + C(1:2)*GL%w(i)*pi
      !d00_inf = d00_mid
      !d00_mid = d00
      !d10_inf = d10_mid
      !d10_mid = d10
      d11_inf = d11_mid
      d11_mid = d11
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels_Lens


subroutine Kernels_Rot(rL,WA,WB,X,kernel)
!Kernels of pol rot reconstruction noise
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(2)
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:)
  !internal
  type(gauss_legendre_params) :: GL
  character(LEN=4) :: selec
  integer :: i, l, lmax
  double precision :: mu, II, al, c1_inv, c2p, c2m, c3
  double precision, dimension(rL(1):rL(2)) :: al0, alp2, alm2
  double precision, dimension(2) :: d00_sup, d00_mid, d00_inf
  double precision, dimension(2) :: ZA22
  double precision, dimension(2) :: ZB22

  !* initialize
  lmax = size(X)
  X = 0d0

  do l = rL(1), rL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp2(l) = dsqrt((al-2)*(al+3))
    alm2(l) = dsqrt((al+2)*(al-1))
  end do

  !* GL quadrature
  call gl_initialize(GL,int((3*max(rL(2),lmax)+1)/2),1d-15)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(kernel)
    case ('Sp','Gp')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(2,2,rL,WB,mu,ZB22)
      II = ZA22(1)*ZB22(1)-ZA22(2)*ZB22(2)
    case ('Sm','Gm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(2,2,rL,WB,mu,ZB22)
      II = ZA22(1)*ZB22(1)+ZA22(2)*ZB22(2)
    case default
      stop 'error: no kernel'
    end select
    do l = 1, lmax ! loop for l
      al = dble(l)
      if(l==1) then 
        d00_sup = wigd_ini(0,0,mu)
        d00_inf = d00_mid
        d00_mid = d00_sup
        d00_sup(1:2) = mu*d00_mid(1:2)*(al*(2d0*al-1d0))/al**2
      else
        c1_inv = (2d0*al-1d0)/al
        c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        d00_sup(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        d00_sup(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
      end if
      X(l) = X(l) + 4d0*II*d00_sup(1)*GL%w(i)*pi
      d00_inf = d00_mid
      d00_mid = d00_sup
    end do
  end do

  !px
  select case(kernel)
  case ('Gp','Gm')
    X = -X
  end select

  call gl_finalize(GL)

end subroutine Kernels_Rot


subroutine Kernels_Tau(rL,WA,WB,X,kernel,gln,gle)
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(2)
  integer, intent(in), optional :: gln
  double precision, intent(in), optional :: gle
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:)
  !internal
  type(gauss_legendre_params) :: GL
  integer :: i, l, lmax, gn
  double precision :: mu, II, al, c1_inv, c3, ge
  double precision, dimension(2) :: d00_sup, d00_mid, d00_inf
  double precision, dimension(2) :: ZA00, ZB00, ZA20, ZB20, ZA22, ZB22

  !* initialize
  lmax = size(X)

  X = 0d0

  !* GL quadrature
  gn = int((3*max(rL(2),lmax)+1)/2)
  ge = 1d-15
  if (present(gln)) gn = gln
  if (present(gle)) ge = gle
  call gl_initialize(GL,gn,ge)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(kernel)
    case ('S0','G0')
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(0,0,rL,WB,mu,ZB00)
      II = 2d0*ZA00(1)*ZB00(1)
    case ('Sc','Gc')
      call ZETA(2,0,rL,WA,mu,ZA20)
      call ZETA(2,0,rL,WB,mu,ZB20)
      II = ZA20(1)*ZB20(1) + ZA20(2)*ZB20(2)
    case ('Sp','Gp')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(2,2,rL,WB,mu,ZB22)
      II = ZA22(1)*ZB22(1) + ZA22(2)*ZB22(2)
    case ('Sm','Gm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(2,2,rL,WB,mu,ZB22)
      II = ZA22(1)*ZB22(1) - ZA22(2)*ZB22(2)
    case default
      stop 'error: no kernel'
    end select
    do l = 1, lmax ! loop for L
      al = dble(l)
      if(l==1) then 
        d00_sup = wigd_ini(0,0,mu)
        d00_inf = d00_mid
        d00_mid = d00_sup
        d00_sup(1:2) = mu*d00_mid(1:2)*(al*(2d0*al-1d0))/al**2
      else
        c1_inv = (2d0*al-1d0)/al
        c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        d00_sup(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        d00_sup(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
      end if
      X(l) = X(l) + II*d00_sup(1)*GL%w(i)*pi
      d00_inf = d00_mid
      d00_mid = d00_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels_Tau


subroutine Kernels_LensTau(rL,WA,WB,X,kernel)
!Kernels of lensing reconstruction x tau
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(2)
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:,:)
  !internal
  type(gauss_legendre_params) :: GL
  integer :: i, l, lmax
  double precision :: pm, mu, II, al, c1_inv, c2p, c2m, c3
  double precision, dimension(rL(1):rL(2)) :: al0, alp2, alm2
  double precision, dimension(2) :: d10, d10_mid, d10_inf
  double precision, dimension(2) :: ZA00, ZA22
  double precision, dimension(2) :: ZB10, ZB32, ZB21

  !* initialize
  lmax = size(X,dim=2)

  if (kernel=='Sp'.or.kernel=='Gp') pm = 1d0
  if (kernel=='Sm'.or.kernel=='Gm') pm = -1d0

  X = 0d0

  do l = rL(1), rL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp2(l) = dsqrt((al-2)*(al+3))
    alm2(l) = dsqrt((al+2)*(al-1))
  end do

  !* GL quadrature
  call gl_initialize(GL,int((3*max(rL(2),lmax)+1)/2),1d-15)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(kernel)
    case ('S0','G0')
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(1,0,rL,WB*al0,mu,ZB10)
      II = ZA00(1)*ZB10(1)*2d0
    case ('Sc','Gc')
    case ('Sp','Gp','Sm','Gm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(3,2,rL,WB*alp2,mu,ZB32)
      call Zeta(2,1,rL,WB*alm2,mu,ZB21)
      II = ((ZB32(1)-ZB21(1))*ZA22(1)+pm*(ZB32(2)+ZB21(2))*ZA22(2))/2d0
    case default
      stop 'error: no kernel'
    end select
    do l = 1, lmax ! loop for l
      al = dble(l)
      if (l==1) then 
        d10_mid  = 0d0 !d^0_10(mu)
        d10 = -dsqrt(1-mu**2)/dsqrt(2.d0) !d^1_10(mu)
      else
        c1_inv = al*(2d0*al-1d0)/(dsqrt((al**2-1d0))*al)
        c3 = dsqrt(((al-1d0)**2-1d0))*(al-1d0)/((al-1d0)*(2d0*al-1d0))
        d10(1) = (mu*d10_mid(1) - c3*d10_inf(1))*c1_inv
        d10(2) = (mu*d10_mid(2) - c3*d10_inf(2))*c1_inv
      end if
      X(1:2,l) = X(1:2,l) + II*d10(1)*dsqrt(l*(l+1d0))*GL%w(i)*pi
      d10_inf = d10_mid
      d10_mid = d10
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels_LensTau


end module alkernel

