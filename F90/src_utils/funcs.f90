!////////////////////////////////////////////////////!
! * Intrinsic Special Functions
!////////////////////////////////////////////////////!

module funcs
  implicit none
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0

  private pi

contains


function factorial(n)  result(f)
  implicit none
  integer, intent(in) :: n
  integer :: i
  double precision :: f

  f = 1d0
  do i = 1, n
    f = f*dble(i)
  end do

end function factorial


!//// Wigner d functions ////!

function wigd_ini(m1,m2,mu) result(d)
!* initial values of the wigner d-function for the recursion formula (ell=min(m1,m2) with m1, m2, mu)
  implicit none
  !I/O
  integer, intent(in) :: m1, m2
  double precision, intent(in) :: mu
  !internal
  double precision :: x(2), si, d(2)

  ! mu = cos(theta)
  x(1) = 1d0 + mu
  x(2) = 1d0 - mu
  si = dsqrt(1d0 - mu**2)

  if (m1==0.and.m2==0) then 
    !(0,0)
    d = 1d0
  else if(m1==1.and.m2==0) then 
    !(1,0)
    d = -si/dsqrt(2.d0)
  else if(m1==1.and.m2==1) then 
    !(1,+/-1)
    d = x*0.5d0
  else if(m1==2.and.m2==0) then 
    !(2,0)
    d = dsqrt(3d0/8d0)*si**2
  else if(m1==2.and.m2==1) then 
    !(2,+/-1)
    d = -si*0.5d0*x
  else if (m1==2.and.m2==2) then 
    !(2,+/-2)
    d = x**2*0.25d0
  else if (m1==3.and.m2==0) then 
    !(3,0)
    d = -dsqrt(5d0)/4d0*si**3
  else if (m1==3.and.m2==1) then 
    !(3,+/-1)
    d(1) = x(1)**2*x(2)*dsqrt(15.d0)*0.125d0
    d(2) = x(2)**2*x(1)*dsqrt(15.d0)*0.125d0
  else if (m1==3.and.m2==2) then 
    !(3,+/-2)
    d = dsqrt(1.5d0)*(-1.d0)*si*x**2/4.d0
  else if (m1==3.and.m2==3) then 
    !(3,+/-3)
    d = x**3*0.125d0
  else
    stop 'error (wigd_ini_factor): no initial values'
  end if

end function wigd_ini


subroutine wigd_recursion(L,M1,M2,mu,d_sup,d_mid,d_inf)
!* recursion formula of the wigner d function
  implicit none
  !I/O
  integer, intent(in) :: L, M1, M2
  double precision, intent(in) :: d_inf(2), d_mid(2), mu
  double precision, intent(out) :: d_sup(2)
  !internal
  double precision :: aL, aM1, aM2, c1_inv, c2p, c2m, c3

  aL = dble(L)
  aM1 = dble(M1)
  aM2 = dble(M2)

  c1_inv = (aL*(2d0*aL-1d0))/dsqrt((aL**2-aM1**2)*(aL**2-aM2**2))

  if(L<=1) then 
    c3 = 0d0
    c2p = - mu
    c2m = - mu
  else
    c2p = aM1*aM2/((aL-1d0)*aL) - mu
    c2m = -aM1*aM2/((aL-1d0)*aL) - mu
    c3 = dsqrt(((aL-1d0)**2-aM1**2)*((aL-1d0)**2-aM2**2))/((aL-1d0)*(2d0*aL-1d0))
  end if

  if (L+1>=M1.and.L+1>=M2) then 
    d_sup(1) = -(c2p*d_mid(1) + c3*d_inf(1))*c1_inv
    d_sup(2) = -(c2m*d_mid(2) + c3*d_inf(2))*c1_inv
  else 
    d_sup = 0d0
  end if

end subroutine wigd_recursion


!//// wigner-3j symbols ////!

function w00_ini(L,L1) result(W00)
!return W3j(L,L1,L1;0,0,0)
  implicit none
  !I/O
  integer, intent(in) :: L, L1
  !internal
  integer :: i
  double precision :: p, L1L0, W00

  !* compute A(L1,L,0) (useful for first values of W_l2)
  L1L0 = 1d0
  do i = 1, L
    L1L0 = L1L0*(dble(i)-0.5d0)*dble(L1+i)/(dble(i)*(dble(L1+i)-0.5d0))
  end do

  !* Wigner 3j (0,0)
  W00 = dsqrt(L1L0)/dsqrt(2d0*dble(L+L1)+1d0)

  !* parity test (L0+L1+L2)
  p = dble(L+L1)*0.5d0
  if (dint(p).ne.p) W00 = -W00 !odd parity

end function w00_ini


function w3j_ini_ratio(L,L1,M1,M2)
!return ratio of W3j(L,L1,L1,M1,M2,-M1-M2)/W3j(L,L1,L1,0,0,0)
  implicit none
  !I/O
  integer, intent(in) :: L, L1, M1, M2
  !internal
  integer :: M3, k1, k2, k3, I1, I2
  double precision :: fac, p
  double precision :: w3j_ini_ratio

  I1 = M1
  I2 = M2
  M3 = -M1-M2

  if(M1<0) I1 = -I1
  if(M2<0) I2 = -I2
  if(M3<0) M3 = -M3

  fac = 1d0
  p = dble(M3)/2d0
  if(dint(p).ne.p)  fac = -1d0

  do k1 = 1, I1
    fac = fac*dsqrt(dble(L1-k1+1)/dble(L1+k1))
  end do
  do k2 = 1, I2
    fac = fac*dsqrt(dble(L-k2+1)/dble(L+k2))
  end do
  do k3 = 1, M3
    fac = fac*dsqrt(dble(L+L1+k3)/dble(L+L1-k3+1))
  end do

  w3j_ini_ratio = fac

end function w3j_ini_ratio


function w3j_sec(L,L1,M1,M2)
  implicit none
  integer, intent(in) :: M1, M2, L, L1
  double precision :: w3j_sec

  w3j_sec = dble(L*M1-L1*M2)*dsqrt(dble(2*L+2*L1+1)/(L*L1*((L+L1)**2-(dble(M1+M2))**2)))

end function w3j_sec


subroutine w3j_recursion(iL0,iL1,iL2,iM1,iM2,W3j)
  implicit none
  !I/O
  integer, intent(in) :: iL0, iL1, iL2, iM1, iM2
  double precision, intent(inout) :: W3j(1:3)
  !internal
  double precision :: L0, L1, L2, M1, M2, a, b1, b2, c

  L0 = dble(iL0)
  L1 = dble(iL1)
  L2 = dble(iL2)
  M1 = dble(iM1)
  M2 = dble(iM2)

  a  = (L2+2)*dsqrt((-L2+L0+L1)*(L2-L0+L1+1)*(L2+L0-L1+1)*(L2+L0+L1+2))
  b1 = 2*(L2+1)*(L2+2)*(2*L2+3)
  b2 = (2*L2+3)*((L2+1)*(L2+2)+L0*(L0+1)-L1*(L1+1))
  c  = (L2+1)*dsqrt((-L2+L0+L1-1)*(L2-L0+L1+2)*(L2+L0-L1+2)*(L2+L0+L1+3))

  if (L2>=M1+M2) then
    W3j(3) = w3j_b(M1,M2,L2,a,b1,b2)*W3j(2) + w3j_c(M1,M2,L2,a,c)*W3j(1)
  else 
    W3j(3) = 0d0
  end if

end subroutine w3j_recursion


function w3j_b(M1,M2,L,a,b1,b2) result(b)
  implicit none
  !I/O
  double precision, intent(in) :: M1, M2, L, a, b1, b2
  !internal
  double precision :: b

  b = -(M2*b1-(M1+M2)*b2)/(dsqrt((L+1)**2-(M1+M2)**2)*a)

end function w3j_b

function w3j_c(M1,M2,L,a,c)
  implicit none
  !I/O
  double precision, intent(in) :: M1, M2, L, a, c
  !internal
  double precision :: w3j_c

  w3j_c = -(dsqrt((L+2)**2-(M1+M2)**2)*c)/(dsqrt((L+1)**2-(M1+M2)**2)*a)

end function w3j_c


!////////////////////////////////////!
! Probability Distribution Functions !

subroutine GET_WISHART_DISTRIBUTION(DOF,n,a1,a2,mean,x,pdf,f)
!output: pdf(n)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: DOF, n
  double precision, intent(in) :: a1, a2, mean
  double precision, intent(out) :: pdf(1:n), x(1:n)
  !intenral
  integer :: i
  double precision :: norm, width, logpdf(n)

  norm  = 0d0
  width = (a2-a1)/(n*mean)
  do i = 1, n
    x(i) = dble(i)*width + a1/mean
    logpdf(i) = - x(i)*dble(DOF) + dble(DOF-1)*log(x(i))
    norm = norm + dexp(logpdf(i))
  end do
  pdf = dexp(logpdf)/norm

  if(present(f)) then 
    open(unit=20,file=trim(f),status='replace')
    do i = 1, n
      write(20,'(3(E12.5,1X))') x(i), pdf(i), logpdf(i)
    end do
    close(20)
  end if

end subroutine GET_WISHART_DISTRIBUTION


!///////////////////!
! Special Functions
! - Hermite polynomials
! - Error Function
! - Gamma Function

function Her(k,x) !Hermite polynomials
  implicit none
  !I/O
  integer, intent(in) :: k
  double precision, intent(in) :: x
  !internal
  integer :: i
  double precision :: Her, y0, y1, y2

  Her = 0d0
  if(k==-1) Her = dsqrt(pi*0.5d0)*dexp(x**2*0.5d0)*erfc(x/dsqrt(2.d0))
  y0 = 1d0
  if(k==0) Her = y0
  y1 = x
  if(k==1) Her = y1
  if(k>=2) then 
    do i = 2, k
      y2 = x*y1-dble(i-1)*y0
      y0 = y1
      y1 = y2
    end do 
    Her = y1
  end if

end function Her


function erfc(x)
  implicit none
  double precision, intent(in) :: x
  double precision erfc

  if (x < 0) then
    erfc = 1.0+GammaFuncP(0.5d0,x*x)
  else
    erfc = GammaFuncQ(0.5d0,x*x)
  end if

end function erfc


function GammaFuncP(a,x)
  implicit none
  double precision, intent(in) :: a,x
  double precision gamser, gammcf, gln, GammaFuncP

  if (x < (a+1.d0)) then
    call gser(gamser,a,x,gln)
    GammaFuncP = gamser
  else
    call gcf(gammcf,a,x,gln)
    GammafuncP = 1.d0 - gammcf
  end if

end function GammaFuncP


function GammaFuncQ(a,x)
  implicit none
  double precision, intent(in) :: a,x 
  double precision gamser, gammcf, gln, GammaFuncQ

  if (x < (a+1d0)) then
    call gser(gamser,a,x,gln)
    GammaFuncQ = 1d0-gamser
  else
    call gcf(gammcf,a,x,gln)
    GammafuncQ = gammcf
  end if

end function GammaFuncQ


subroutine gser(gamser,a,x,gln)
  implicit none
  double precision, intent(in) :: a,x
  double precision, intent(out) :: gamser, gln
  double precision :: EPS, ap, sum, del
  integer :: ITMAX, n

  ITMAX = 100
  EPS   = 3d-7 !numerical error 
  gln   = lnGamma(a)

  if (x<=0d0) then
    if (x < 0) write(*,*) 'x less than 0 in routine gser'
    gamser = 0d0
  else
    ap = a
    sum = 1d0/a
    del = sum
    do n = 1, ITMAX
      ap = ap + 1d0
      del = del*(x/ap)
      sum = sum + del
      if (abs(del)<abs(sum)*EPS) gamser = sum*dexp(-x+a*dlog(x)-gln)
    end do
  end if

end subroutine gser


subroutine gcf(gammcf,a,x,gln)
  implicit none
  double precision, intent(in) :: a,x
  double precision, intent(out) :: gammcf, gln
  double precision :: EPS, FPMIN, gold, fac, g, b1, b0, anf, ana, an, a1, a0
  integer :: ITMAX, n

  ITMAX = 100
  EPS   = 3d-7 !numerical error
  b0    = 0d0
  gold  = 0d0
  fac   = 1d0
  b1    = 1d0
  a0    = 1d0
  gln   = lnGamma(a)

  a1 = x
  do n = 1, ITMAX
    an = dble(n)
    ana = an -a
    a0 = (a1+a0*ana)*fac
    b0 = (b1+b0*ana)*fac
    anf = an*fac
    a1 = x*a0+anf*a1
    b1 = x*b0+anf*b1
    if (.not.a1==0) then
      fac = 1d0/a1
      g = b1*fac
      if (abs(g-gold)/g < EPS) then
        gammcf = dexp(-x+a*dlog(x)-gln)*g
        exit
      end if
      gold = g
    end if
  end do

end subroutine gcf


function lnGamma(x) !Lanczos formula (gamma=5, N=6)
  implicit none
  double precision, intent(in) :: x
  double precision tmp, ser, lnGamma
  double precision :: cof(6), cof_0
  integer j

  cof_0  = 1.000000000190015d0
  cof(1) = 76.18009172947146d0
  cof(2) = -86.50532032941677d0
  cof(3) = 24.01409824083091d0
  cof(4) = -1.231739572450155d0
  cof(5) = 0.1208650973866179d-2
  cof(6) = -0.5395239384953d-5

  tmp = (x+0.5d0)*dlog(x+5.5d0) - (x+5.5d0)
  ser = cof_0
  do j = 1, 6
    ser = ser + cof(j)/(x+j)
  end do
  lnGamma = tmp + dlog(2.5066282746310005d0*ser/x) !ln[Gammma(x+1)]-ln(x)

end function lnGamma


!//// Bessel functions ////!

function Bessel_J(n,x,eps)  result(f)
! Bessel function
  implicit none
  !I/O
  integer, intent(in) :: n
  integer, intent(in), optional :: eps
  double precision, intent(in) :: x
  !internal
  integer :: i, m
  double precision :: f, dt, t

  m = 32
  if(present(eps)) m = eps

  t = 0d0
  dt = 2d0*pi/dble(m)
  f = 0d0
  do i = 1, m-1
    t = t + dt
    f = f + dcos(x*dsin(t)-n*t)*dt/(2d0*pi)
  end do
  f = f + 1d0*dt/(2d0*pi)

end function Bessel_J



!//// HyperGeometric Function ////!

subroutine hygfx(a,b,c,x,hf)
! //// Hypergeometric function F(A,B,C,X) //// !
! * Licensing:
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.
!
! * Modified:
!    08 September 2007
!
  implicit none
  integer(kind=4) :: j,k,n,nm,m
  logical :: l0,l1,l2,l3,l4,l5
  double precision :: a,a0,aa,b,bb,c,c0,c1,eps,f0,f1,g0,g1,g2,g3,ga,gabc,gam,sp,x1
  double precision :: gb,gbm,gc,gca,gcab,gcb,gm,hf,hw,pa,pb,r,r0,r1,rm,rp,sm,sp0,x
  double precision, parameter :: el = 0.5772156649015329

  l0 = ( c == aint(c) ) .and. ( c < 0d0 )
  l1 = ( 1d0 - x < 1.d-15 ) .and. ( c-a-b <= 0d0 )
  l2 = ( a == aint(a) ) .and. ( a < 0d0 )
  l3 = ( b == aint(b) ) .and. ( b < 0d0 )
  l4 = ( c - a == aint(c-a) ) .and. ( c-a <= 0d0 )
  l5 = ( c - b == aint(c-b) ) .and. ( c-b <= 0d0 )

  if ( l0 .or. l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYGFX - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    return
  end if

  if ( 0.95D0 < x ) then 
    eps = 1.0D-8
  else
    eps = 1.0D-15
  end if

  if ( x == 0d0 .or. a == 0d0 .or. b == 0d0 ) then
    hf = 1d0
    return
  else if ( 1d0-x == eps .and. 0d0 < c-a-b) then
    call gamma(c,gc)
    call gamma(c-a-b,gcab)
    call gamma(c-a,gca)
    call gamma(c-b,gcb)
    hf = gc*gcab/(gca*gcb)
    return
  else if (1d0+x<=eps .and. abs(c-a+b-1d0)<=eps) then
    g0 = sqrt(pi)*2d0**(-a)
    call gamma(c,g1)
    call gamma(1d0 + a/2d0 - b, g2 )
    call gamma(0.5d0 + 0.5d0*a, g3 )
    hf = g0*g1/(g2*g3)
    return
  else if ( l2 .or. l3 ) then
    if(l2) nm = int(abs(a))
    if(l3) nm = int(abs(b))
    hf = 1d0
    r  = 1d0
    do k = 1, nm
      r  = r*(a+k-1d0)*(b+k-1d0)/(k*(c+k-1d0))*x
      hf = hf + r
    end do
    return
  else if ( l4 .or. l5 ) then
    if(l4) nm = int(abs(c-a))
    if(l5) nm = int(abs(c-b))
    hf = 1d0
    r  = 1d0
    do k = 1, nm
      r = r*(c-a+k-1d0)*(c-b+k-1d0)/(k*(c+k-1d0))*x
      hf = hf + r
    end do
    hf = (1d0-x)**(c-a-b)*hf
    return
  end if
  aa = a
  bb = b
  x1 = x
!
!  WARNING: ALTERATION OF INPUT ARGUMENTS A AND B, WHICH MIGHT BE CONSTANTS.
!
  if ( x < 0d0 ) then
    x = x/(x-1d0)
    if ( a < c .and. b < a .and. 0d0 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if
  if ( 0.75D0 <= x ) then
    gm = 0d0
    if ( abs(c-a-b-aint(c-a-b)) < 1d-15 ) then
      m = int(c-a-b)
      call gamma(a,ga)
      call gamma(b,gb)
      call gamma(c,gc)
      call gamma(a+m,gam)
      call gamma(b+m,gbm)
      call psi(a,pa)
      call psi(b,pb)
      if (m/= 0) gm = 1d0
      do j = 1, abs(m)-1
        gm = gm*j
      end do
      rm = 1d0
      do j = 1, abs(m)
        rm = rm * j
      end do
      f0 = 1d0
      r0 = 1d0
      r1 = 1d0
      sp0 = 0d0
      sp = 0d0
      if ( 0 <= m ) then
        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1d0 )**m / ( ga * gb * rm )
        do k = 1, m - 1
          r0 = r0*(a+k-1d0)*(b+k-1d0)/(k*(k-m))*(1d0-x)
          f0 = f0 + r0
        end do
        do k = 1, m
          sp0 = sp0+1d0/(a+k-1d0)+1d0/(b+k-1d0)-1d0/dble(k)
        end do
        f1 = pa + pb + sp0 + 2d0*el + log(1d0-x)
        hw = f1
        do k = 1, 250
          sp = sp + (1d0-a)/(k*(a+k-1d0))+(1d0-b)/(k*(b+k-1d0))
          sm = 0d0
          do j = 1, m
            sm = sm+(1d0-a)/((j+k)*(a+j+k-1d0))+1d0/(b+j+k-1d0)
          end do
          rp = pa + pb + 2d0*el + sp + sm + log(1d0-x)
          r1 = r1*(a+m+k-1d0)*(b+m+k-1d0)/(k*(m+k))*(1d0-x)
          f1 = f1 + r1 * rp
          if ( abs(f1-hw) < abs(f1)*eps )  exit
          hw = f1
        end do
        hf = f0 * c0 + f1 * c1
      else if ( m < 0 ) then
        m = - m
        c0 = gm * gc / ( ga * gb * ( 1d0 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )
        do k = 1, m - 1
          r0 = r0*(a-m+k-1d0)*(b-m+k-1d0)/(k*(k-m))*(1d0-x)
          f0 = f0 + r0
        end do
        do k = 1, m
          sp0 = sp0 + 1d0/dble(k)
        end do
        f1 = pa + pb - sp0 + 2d0*el + log(1d0-x)
        do k = 1, 250
          sp = sp + (1d0-a)/(k*(a+k-1d0)) + (1d0-b)/(k*(b+k-1d0))
          sm = 0d0
          do j = 1, m
            sm = sm + 1d0/dble(j+k)
          end do
          rp = pa+pb+2d0*el+sp-sm+log(1d0-x)
          r1 = r1*(a+k-1d0)*(b+k-1d0)/(k*(m+k))*(1d0-x)
          f1 = f1+r1*rp
          if ( abs(f1-hw) < abs(f1)*eps )  exit
          hw = f1
        end do
        hf = f0 * c0 + f1 * c1
      end if
    else
      call gamma(a,ga)
      call gamma(b,gb)
      call gamma(c,gc)
      call gamma(c-a,gca)
      call gamma(c-b,gcb)
      call gamma(c-a-b,gcab)
      call gamma(a+b-c,gabc)
      c0 = gc*gcab/(gca*gcb)
      c1 = gc*gabc/(ga*gb) * (1d0-x)**(c-a-b)
      hf = 0d0
      r0 = c0
      r1 = c1
      hw = hf !namikawa
      do k = 1, 250
        r0 = r0*(a+k-1d0)*(b+k-1d0)/(k*(a+b-c+k))*(1d0-x)
        r1 = r1*(c-a+k-1d0)*(c-b+k-1d0)/(k*(c-a-b+k))*(1d0-x)
        hf = hf + r0 + r1
        if ( abs(hf-hw) < abs(hf)*eps )   exit
        hw = hf
      end do
      hf = hf + c0 + c1
    end if
  else
    a0 = 1d0
    if (a<c .and. c<2d0*a .and. b<c .and. c<2d0*b ) then
      a0 = (1d0-x)**(c-a-b)
      a = c-a
      b = c-b
    end if
    hf = 1d0
    r  = 1d0
    hw = hf !namikawa
    do k = 1, 250
      r = r*(a+k-1d0)*(b+k-1d0)/(k*(c+k-1d0))*x
      hf = hf+r
      if ( abs(hf-hw) <= abs(hf)*eps )  exit
      hw = hf
    end do
    hf = a0 * hf
  end if
  if ( x1 < 0d0 ) then
    x = x1
    c0 = 1d0 / ( 1d0 - x )**aa
    hf = c0 * hf
  end if
  a = aa
  b = bb
  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYGFX - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if
  return

end subroutine hygfx


subroutine gamma(x,ga)
! //// Gamma function //// !
! * Licensing:
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.
! * Modified:
!    08 September 2007
  implicit none
  double precision, dimension(26) :: g = (/1.0D0, 0.5772156649015329D0, &
   -0.6558780715202538D0, -0.420026350340952D-1, 0.1665386113822915D0, &
   -0.421977345555443D-1, -0.96219715278770D-2, 0.72189432466630D-2, &
   -0.11651675918591D-2, -0.2152416741149D-3, 0.1280502823882D-3, & 
   -0.201348547807D-4, -0.12504934821D-5, 0.11330272320D-5, &
   -0.2056338417D-6, 0.61160950D-8, 0.50020075D-8, -0.11812746D-8, &
    0.1043427D-9, 0.77823D-11, -0.36968D-11, 0.51D-12, &
   -0.206D-13, -0.54D-14, 0.14D-14, 0.1D-15 /)
  integer(kind=4) :: k,m,m1
  double precision :: r,x,z,ga,gr

  if (x == aint(x)) then
    if (0.0D0 < x) then
      ga = 1.0D0
      m1 = int(x) - 1
      do k = 2, m1
        ga = ga*k
      end do
    else
      ga = 1.0D+300
    end if
  else
    if ( 1.0D0 < abs(x) ) then
      z = abs(x)
      m = int(z)
      r = 1.0D0
      do k = 1, m
        r = r*(z-real(k,kind=8))
      end do
      z = z-real(m,kind=8)
    else
      z = x
    end if
    gr = g(26)
    do k = 25, 1, -1
      gr = gr*z + g(k)
    end do
    ga = 1.0D0/(gr*z)
    if ( 1.0D0 < abs(x) ) then
      ga = ga*r
      if ( x < 0.0D0 ) then
        ga = -pi/(x*ga*sin(pi*x))
      end if
    end if
  end if
  return

end subroutine gamma


subroutine psi(x,ps)
! //// PSI function //// !
! * Licensing:
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.
!
!  Modified:
!    08 September 2007
!
  implicit none
  integer(kind=4) :: k, n
  double precision :: ps,s,x,x2,xa
  double precision, parameter :: a1 = -0.083333333333333333
  double precision, parameter :: a2 =  0.0083333333333333333
  double precision, parameter :: a3 = -0.0039682539682539683
  double precision, parameter :: a4 =  0.0041666666666666667
  double precision, parameter :: a5 = -0.0075757575757575758
  double precision, parameter :: a6 =  0.021092796092796093
  double precision, parameter :: a7 = -0.083333333333333333
  double precision, parameter :: a8 =  0.4432598039215686
  double precision, parameter :: el = 0.5772156649015329

  xa = abs(x)
  s = 0.0D0

  if ( x == aint(x) .and. x <= 0.0D0 ) then
    ps = 1.0D+300
    return
  else if ( xa == aint(xa) ) then
    n = int (xa)
    do k = 1, n - 1
      s = s + 1.0D0/real(k,kind=8)
    end do
    ps = - el + s
  else if ( xa + 0.5D0 == aint(xa+0.5D0) ) then
    n = int(xa-0.5D0)
    do k = 1, n
      s = s + 1.0D0 / real(2*k-1,kind=8)
    end do
    ps = - el + 2.0D0 * s - 1.386294361119891D0
  else
    if ( xa < 10.0D0 ) then
      n = 10 - int(xa)
      do k = 0, n - 1
        s = s + 1.0D0 / ( xa + real(k,kind=8) )
      end do
      xa = xa + real(n,kind=8)
    end if
    x2 = 1.0D0/(xa*xa)
    ps = log(xa) - 0.5D0/xa + x2*((((((( a8 * x2 + a7 ) * x2 + a6 ) &
      * x2 + a5 ) * x2 + a4 ) * x2 + a3 ) * x2 + a2 ) * x2 + a1 )
    ps = ps - s
  end if
  if ( x < 0.0D0 ) then
    ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D0 / x
  end if
  return

end subroutine psi


end module funcs

