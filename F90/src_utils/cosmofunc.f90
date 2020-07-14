!////////////////////////////////////////////////////!
! * My functions
!////////////////////////////////////////////////////!

module cosmofunc
  use funcs,   only: hygfx, erfc, Her
  use general, only: gl_initialize, gl_finalize, gauss_legendre_params, GLdxs, GLpoint, neighb
  implicit none

  interface C_z
    module procedure C_z_single, C_z_array
  end interface

  interface D_z
    module procedure D_z_single, D_z_array
  end interface

  interface g_z
    module procedure g_factor, g_factor_array
  end interface

  interface H_z
    module procedure H_z_single, H_z_array
  end interface

  !* cosmological parameters
  type cosmoparams !nu=On/Om
    double precision :: H0=70d0, h=0.7d0, Ov=0.7d0, Om=0.3d0, w0=-1d0, wa=0d0, nu=4.38225537d-2, ns=0.9645d0
  end type cosmoparams

  !* local parameters
  double precision, parameter :: c = 2.99792458d5
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0
  private c, pi

  private hygfx, erfc, Her
  private gl_initialize, gl_finalize, gauss_legendre_params, GLdxs, GLpoint, neighb

contains

subroutine set_cosmoparams(cp,cpmodel)
  ! cosmological parameters
  implicit none
  character(*), intent(in) :: cpmodel
  type(cosmoparams), intent(inout) :: cp

  select case(cpmodel)
  case ('model0')
    cp%Om = 0.31201550908d0
    cp%H0 = 67.51d0
    cp%w0 = -1d0
    cp%wa = 0d0
    cp%nu = 0.06d0/(93.14d0*(cp%H0/100d0)**2*cp%Om)
  case ('modelw')
    cp%Om = 0.279d0
    cp%H0 = 70d0
    cp%w0 = -1d0
    cp%wa = 0d0
    cp%nu = 0d0
  case ('modelp')
    cp%Om = 0.3156d0
    cp%H0 = 67.27d0
    cp%w0 = -1d0
    cp%wa = 0d0
    cp%nu = 0.00443d0
  end select
  cp%h  = cp%H0/100d0
  cp%Ov = 1d0 - cp%Om
  cp%ns = 0.9645d0

end subroutine set_cosmoparams

! for aps calculation
function kernel_cgg(z,cp)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z(:)
  integer :: i
  double precision :: f(size(z))

  do i = 1, size(z)
    f(i) = H_z(z(i),cp) / (C_z(z(i),cp))**2
  end do

end function kernel_cgg


function omega_m(a,cp) 
! Omega_m(a) = rho_m(a) / rho_c(a)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: a
  double precision :: omega_m, omega_t, Qa2, Ok

  Ok = 1d0-cp%Om-cp%Ov
  Qa2     = a**(-1d0-3d0*(cp%w0+cp%wa))*dexp(-3d0*(1d0-a)*cp%wa)
  omega_t = 1d0 - Ok/(Ok+cp%Ov*Qa2+cp%Om/a)
  omega_m = omega_t*cp%Om/(cp%Om+cp%Ov*a*Qa2)

end function omega_m


function omega_v(a,cp) 
! Omega_v = rho_v(a) / rho_c(a)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: a
  double precision :: omega_v, omega_t, Qa2, Ok

  Ok = 1d0-cp%Om-cp%Ov
  Qa2     = a**(-1d0-3d0*(cp%w0+cp%wa))*dexp(-3d0*(1d0-a)*cp%wa)
  omega_t = 1d0 - Ok/(Ok+cp%Ov*Qa2+cp%Om/a)
  omega_v = omega_t*cp%Ov*Qa2/(cp%Ov*Qa2+cp%Om/a)

end function omega_v


function wpca(a,wi,dz,b)  result(w)
  implicit none
  double precision, intent(in) :: a, wi(:), dz, b
  integer :: i
  double precision :: w, z, wfunc

  w = -1d0
  z = 1d0/a -1d0
  do i = 1, size(wi)
    wfunc = atan((z-dz*dble(i-1))*b) - atan((z-dz*dble(i))*b)/pi
    w = w + (wi(i)+1d0)*wfunc
  end do

end function wpca


function rho_de(a,cp,wpca,dz)  result(f)
! * energy density of DE
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: a
  double precision, intent(in), optional :: wpca(:), dz
  integer :: i , N
  double precision :: f, z

  z = 1d0/a - 1d0
  if (present(wpca).and.present(dz)) then
    do i = 1, size(wpca)
      if ( dz*dble((i-1)) <= z .and. z< dz*dble(i)) N = i
    end do
    if(N>1) then
      f = cp%Ov
      do i = 1, N-1
        f = f*(1d0+dz*dble(i))**(3d0*(wpca(i)-wpca(i+1)))
      end do
      f = f*a**(-3d0*(1d0+wpca(N)))
    else !outer region
      f = f*a**(-3d0*(1d0+cp%w0))
    end if
  else
    f = cp%Ov*a**(-3d0*(1d0+cp%w0+cp%wa))*dexp(-3d0*cp%wa*(1d0-a))
  end if

end function rho_de


function drho_da_de(a,cp)  result(f)  
! * d(rho)/da
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: a
  integer :: i , N
  double precision :: f

  f = rho_de(a,cp)*((-3d0*(1d0+cp%w0+cp%wa))/a+3d0*cp%wa)

end function drho_da_de


function g_factor(z,Om,Ov,wz)
! Growth factor at linear perturbation. 
! g(a) = a in the Einstein de Sitter universe. 
! See Eq.(3.8) of arXiv:1105.4825 for the algorithm
  implicit none
!
! [inputs]
! z  --- redshift
! Om --- Omega_m at z=0
! Ov --- Omega_v at z=0
! wz --- w(a) at input z
  double precision, intent(in) :: z, Om, Ov, wz
!
! [internal]
  double precision :: g_factor, a, b, c, zz, gz

  a = -1d0 / (3d0 * wz)
  b = (wz - 1d0)/ (2d0 * wz)
  c = 1d0 - 5d0 / (6d0 * wz)
  zz = (-Ov/Om) * (1d0+z)**(3d0*wz)
  call HYGFX(a,b,c,zz,gz)
  g_factor = gz/(1d0+z)

end function g_factor


function g_factor_array(z,cp,n)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  integer, intent(in) :: n
  double precision, intent(in) :: z(n)
  integer :: i
  double precision :: f(n), w

  f = 1d0
  do i = 1, n
    if (z(i)==0d0) cycle
    w = cp%w0 + (1d0-1d0/(1d0+z(i)))*cp%wa
    f(i) = g_factor(z(i),cp%Om,cp%Ov,w)
  end do

end function g_factor_array


function g_rate(z,cp)
! Growth rate at linear perturbation. 
! f(a) = 1 in the Einstein de Sitter universe. 
! See Eq.(A.9) of arXiv:1105.4825 for the algorithm
  implicit none
!
! [inputs]
! z  --- redshift
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z
!
! [internal]
  double precision :: g_rate
  double precision :: wz, a, b, c, zz, gz0, gz1, Omp

  !g_rate = dlog(g_factor(z,Om,Ov,wz))/dlog(a)

  wz = cp%w0 + (z/(1d0+z))*cp%wa 

  a = -1d0 / (3d0 * wz)
  b = (wz - 1d0)/ (2d0 * wz)
  c = 1d0 - 5d0 / (6d0 * wz)
  zz = (-cp%Ov/cp%Om) * (1d0+z)**(3d0*wz)
  call HYGFX(a,b,c,zz,gz0)
  call HYGFX(a+1,b+1,c+1,zz,gz1)

  ! If F(a) = a^-f(a), dF(a)/da = - F(a) [ f(a)/a + ln(a)*df/da ]
  ! Here, f(a) = 3w(a)
  Omp = -zz * ( 3d0*wz*(1d0+z) + 3d0*dlog(1d0+z)*cp%wa ) ! dOm/da / Om^2 = du/da
  g_rate = 1 + a*b/(1+z)/c * Omp * gz1/gz0

end function g_rate


function D_z_single(z,cp)  result(f)
! Linear growth factor, normalized to the current value
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z
  double precision :: a, w, f

  a = 1d0/(1d0+z)
  w = cp%w0 + (1d0-a)*cp%wa
  f = g_factor(z,cp%Om,cp%Ov,w)/g_factor(0d0,cp%Om,cp%Ov,w)

end function D_z_single


function D_z_array(z,cp,n)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  integer, intent(in) :: n
  double precision, intent(in) :: z(n)
  integer :: i
  double precision :: f(n)

  f = 1d0
  do i = 1, n
    if (z(i)==0d0) cycle
    f(i) = D_z_single(z(i),cp)
  end do

end function D_z_array


!//// Cosmological Distance ////!

function C_z_single(z,cp,n)  result(f)
! * computing the comoving distance
! 
  implicit none
! [inputs]
!   z  --- redshift
!   cp --- containing [H0,OL,Om,w0,wa]
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z
!
! (optional)
!   n  --- number of GL points
  integer, intent(in), optional :: n
!
! [internal]
  type(gauss_legendre_params) :: GL
  integer :: i, gln
  double precision :: f, zi, dz

  gln = int(z+5d0)
  if (present(n)) gln = n
  call gl_initialize(GL,gln,eps=1d-6)

  f = 0d0
  do i = 1, GL%n
    dz = 0.5d0*z*GL%w(i)
    zi = 0.5d0*z*(1d0+GL%z(i))
    f = f + dz/H_z(zi,cp)
  end do

  call gl_finalize(GL)

end function C_z_single


function C_z_array(z,cp)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z(:)
  integer :: i
  double precision :: f(size(z))

  f = 0d0
  do i = 1, size(z)
    f(i) = C_z_single(z(i),cp)
  end do

end function C_z_array


function r_s(z,Omh,Obh,Orh,n)  result(f)
! * computing the sound horizon
! 
  implicit none
! [inputs]
!   z   --- redshift
!   Omh --- Omega_mh^2
!   Obh --- Omega_bh^2
!   Orh --- Omega_rh^2
  double precision, intent(in) :: z, Omh, Obh, Orh
!
! (optional)
!   n       --- number of GL points
  integer, intent(in), optional :: n
!
! [internal]
  type(gauss_legendre_params) :: GL
  integer :: i, gln
  double precision :: f, R, aeq, cs, a, a_s, da

  gln = int(z+5d0)
  if (present(n)) gln = n
  call gl_initialize(GL,gln,eps=1d-6)

  R   = Obh/Orh
  aeq = Orh/Omh
  a_s = 1d0/(1d0+z)
  f = 0d0
  do i = 1, GL%n
    da = 0.5d0*a_s*GL%w(i)
    a  = 0.5d0*a_s*(1d0+GL%z(i))
    cs = 1d0/dsqrt(3d0*(1d0+a*R))
    f = f + c*cs*da/dsqrt(Omh*1d4*(a+aeq))
  end do

  call gl_finalize(GL)

end function r_s


function L_z(z,Cz,cp)  result(f)
! * Luminosity distance
  implicit none
  type(cosmoparams), intent(in), optional :: cp
  double precision, intent(in) :: z
  double precision, intent(in), optional :: Cz
  double precision :: f

  if (present(cp)) f = C_z(z,cp)*(1d0+z)
  if (present(Cz)) f = Cz*(1d0+z)

end function L_z


function dL_dz(z,CzHz,cp)  result(f)
  implicit none
  type(cosmoparams), intent(in), optional :: cp
  double precision, intent(in) :: z
  double precision, optional :: CzHz(1:2)
  double precision :: f

  if (present(cp))   f = C_z(z,cp) + (1d0+z)/H_z(z,cp)
  if (present(CzHz)) f = CzHz(1) + (1d0+z)/CzHz(2)

end function dL_dz


function H_z_single(z,cp)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z
  double precision :: f, a, Ok

  a  = 1d0/(1d0+z)
  Ok = 1d0-cp%Om-cp%Ov
  if (cp%w0==-1d0.and.cp%wa==0d0) then
    f = dsqrt(Ok/a**2+cp%Om/a**3+cp%Ov)
  else
    f = dsqrt(Ok/a**2+cp%Om/a**3+rho_de(a,cp))
  end if
  f = (cp%H0/c)*f

end function H_z_single


function H_z_array(z,cp)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z(:)
  integer :: i
  double precision :: f(size(z))

  f = 0d0
  do i = 1, size(z)
    f(i) = H_z_single(z(i),cp)
  end do

end function H_z_array


function dH_dz(z,cp)  result(f)
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: z
  double precision :: f, a, Ok

  Ok = 1d0 - cp%Om - cp%Ov
  a  = 1d0/(1d0+z)
  if (cp%w0==-1d0.and.cp%wa==0d0) then
    f = (cp%H0/c)**2 * ( Ok*(1d0+z) + 1.5d0*cp%Om*(1d0+z)**2 ) / H_z(z,cp)
  else
    f = (cp%H0/c)**2 * ( Ok*(1d0+z) + 1.5d0*cp%Om*(1d0+z)**2 - 0.5d0*(1d0+z)**2*drho_da_de(a,cp) ) / H_z(z,cp)
  end if

end function dH_dz


function zoftau(tau,tau0,cp,eps) result(f)
! * Return z from tau (z<10)
! * This subroutine utilize the Newton method
! * The initial value is given as an approximate formula
!    zi ~ dt/(dh-dt/3)
!   where dt=tau0-tau and dh=c/H0. 
! * Typycally, the do loop is finished witin 1 (2) times at z<1 (z>1). 
!
  implicit none
! [inputs]
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: tau, tau0
!
! (optional)
  double precision, intent(in), optional :: eps
!
! [internal]
  integer :: i, n, itr
  double precision :: f, fz, pz, ac, Dt, dh
  double precision, allocatable :: t(:), z(:)

  Dt = tau0 - tau
  dh = c/cp%H0
  n  = 10
  ac = 1d-3
  if (present(eps)) ac = eps

  allocate(t(n),z(n));  t=0d0;  z=0d0
  z(1) = Dt/(dh-Dt/3d0)  ! a solution assuming C_z ~ dh*z/(1+z/3)
  f = z(1)
  itr = 0
  if (z(1)>1d-2) then
    ! Newton method
    do i = 1, n-1
      fz = C_z(z(i),cp) - Dt
      pz = 1d0/H_z(z(i),cp)
      z(i+1) = z(i) - fz/pz
      if(z(i+1)<0d0) z(i+1)=0d0
      f  = z(i+1) 
      if (abs(z(i+1)/z(i)-1d0)<ac) exit
      itr = itr + 1
    end do
  end if
  if (itr==n) stop "error: z(tau) does not converge"
  deallocate(t,z)

end function zoftau


!//// Transfer function and power spectrum ////!

function bbksfit(x)  result(f)
! BBKS fitting formula
  implicit none
  double precision, intent(in) :: x
  double precision :: s, f

  s = 1.58113883d0 !sqrt(2.5)
  f = 0.5d0*x*(x**2-3d0)*(erfc(s*x)+erfc(s*x*0.5d0)) &
      + dsqrt(0.4d0/pi)*((7.45d0*x**2+1.6d0)*dexp(-x**2/1.6d0) &
      + (x**2*0.5d0-1.6d0)*dexp(-2.5d0*x**2))

end function bbksfit


function EisensteinHu(k,Tcmb,h,omegab,omegac)  result(f)
! No-wiggle P(k)
  implicit none
  double precision, intent(in) :: k, Tcmb, h, omegab, omegac
  double precision :: f, r, q, Oh, Ob, L0, C0, g, a, s

  Ob = omegab*h**2
  Oh = (omegab + omegac)*h**2
  r = Ob/Oh
  s = 44.5d0*dlog(9.83d0/Oh)/sqrt(1d0+10d0*Oh**(0.75d0))
  a = 1d0 - 0.328d0*dlog(431d0*Oh)*r + 0.38d0*dlog(22.3d0*Oh)*r**2
  g = (Oh/h)*(a + (1d0-a)/(1d0+(0.43d0*k*s*h)**4))
  q = k*(Tcmb/2.7d0)**2/g
  C0 = 14.2d0 + 731d0/(1d0+62.5d0*q)
  L0 = dlog(2d0*dexp(1d0)+1.8d0*q)

  f = L0/(L0 + C0*q**2)

end function EisensteinHu


function neffective(k,ns,Tcmb,h,omegab,omegac) result(n)
  implicit none
  double precision, intent(in) :: k, ns, Tcmb,h,omegab,omegac
  double precision :: n, dk, kp, km, Tp, Tm

  dk = k*1d-4
  kp = k + dk
  km = k - dk
  Tp = EisensteinHu(kp,Tcmb,h,omegab,omegac)
  Tm = EisensteinHu(km,Tcmb,h,omegab,omegac)
  n = 2d0*dlog(Tp/Tm)/dlog(kp/km) + ns

end function neffective


subroutine Peacock_Dodds(rk,rn,rknl,plin,pnl,g)
  implicit none
  double precision, intent(in) :: rk, rn, plin, g
  double precision, intent(out) :: rknl, pnl
  double precision :: A, B, alpha, beta, V, f, n, y

  y = plin
  n = 1d0 + rn/3d0
  A = 0.482d0*n**(-0.947)
  B = 0.226d0*n**(-1.778)
  alpha = 3.31d0*n**(-0.244)
  beta = 0.862d0*n**(-0.287)
  V = 11.55d0*n**(-0.423)

  f = y*((1d0+B*beta*y+(A*y)**(alpha*beta))/(1d0+((A*y)**alpha*g**3/V/dsqrt(y))**beta))**(1.d0/beta)

  pnl = f
  rknl = (1d0+f)**(1d0/3d0)*rk

end subroutine Peacock_Dodds


subroutine halofit_S02(rk,rn,rncur,rknl,plin,pnl,Om,Ov)
! Halo model nonlinear fitting formula as described in Appendix C of Smith et al. (2002)
  implicit none
  double precision, intent(in) :: rk,rn,rncur,rknl,plin,Om,Ov
  double precision, intent(out) :: pnl
  double precision :: gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
  double precision :: y,f1a,f2a,f3a,f1b,f2b,f3b,frac,pq,ph

  gam = 0.86485+0.2989*rn+0.1631*rncur
  a = 1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+0.1670756*rn*rn*rn*rn-0.620695*rncur
  a = 10**a
  b = 10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
  c = 10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
  xmu = 10**(-3.54419+0.19086*rn)
  xnu = 10**(0.95897+1.2857*rn)
  alpha = 1.38848+0.3701*rn-0.1452*rn*rn
  beta = 0.8291+0.9854*rn+0.3400*rn**2
  if(abs(1.d0-Om).gt.0.01) then ! omega evolution 
    f1a = Om**(-0.0732)
    f2a = Om**(-0.1423)
    f3a = Om**(0.0725)
    f1b = Om**(-0.0307)
    f2b = Om**(-0.0585)
    f3b = Om**(0.0743)       
    frac = Ov/(1.-Om) 
    f1 = frac*f1b + (1-frac)*f1a
    f2 = frac*f2b + (1-frac)*f2a
    f3 = frac*f3b + (1-frac)*f3a
  else         
    f1 = 1.0
    f2 = 1.
    f3 = 1.
  end if
  y = (rk/rknl)
  ph = a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
  ph = ph/(1+xmu*y**(-1)+xnu*y**(-2))
  pq = plin*(1+plin)**beta/(1+plin*alpha)*exp(-y/4.0-y**2/8.0)
  pnl = pq+ph

end subroutine halofit_S02


subroutine halofit_T12(rk,rn,rncur,rknl,plin,pnl,Om_m,Om_v,Omm0,w,fnu)
! - Taken from CAMB
! - plin is the dimenssion less power spectrum
!
  implicit none
  double precision, intent(in) :: rk,rn,rncur,rknl,plin,Om_m,Om_v,Omm0,w,fnu
  double precision, intent(out) :: pnl
  double precision :: gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
  double precision :: pq,ph,plinaa
  double precision :: y,f1a,f2a,f3a,f1b,f2b,f3b,frac

  gam = 0.1971-0.0843*rn+0.8460*rncur
  a = 1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*Om_v*(1.+w)
  a = 10**a
  b = 10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+0.2279*om_v*(1.+w))
  c = 10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
  xmu=0.
  xnu=10**(5.2105+3.6902*rn)
  alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
  beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+0.3980*rn**4-0.1682*rncur + fnu*(1.081 + 0.395*rn**2)

  if(abs(1.d0-Om_m).gt.0.01) then ! omega evolution 
    f1a=Om_m**(-0.0732)
    f2a=Om_m**(-0.1423)
    f3a=Om_m**(0.0725)
    f1b=Om_m**(-0.0307)
    f2b=Om_m**(-0.0585)
    f3b=Om_m**(0.0743)       
    frac=Om_v/(1.-Om_m) 
    f1=frac*f1b + (1-frac)*f1a
    f2=frac*f2b + (1-frac)*f2a
    f3=frac*f3b + (1-frac)*f3a
  else         
    f1=1.0
    f2=1.
    f3=1.
  endif
  y=(rk/rknl)

  ph = a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
  ph = ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+fnu*0.977)
  plinaa = plin*(1+fnu*47.48*rk**2/(1+1.5*rk**2))
  pq = plin*(1+plinaa)**beta/(1+plinaa*alpha)*exp(-y/4.0-y**2/8.0)
  pnl = pq + ph

end subroutine halofit_T12


subroutine halofit_nu(k,CC,CB,BB,NN,CN,BN,n,nc,knl,fc,fb,fcb,Pnl,Om,Ov)
! Halofit with linear neutrinos
  implicit none
  double precision, intent(in) :: n, nc, knl, fc, fb, fcb, Om, Ov
  double precision, intent(in) :: k, CC, CB, BB, NN, CN, BN
  double precision, intent(out) :: pnl
  double precision :: Plin, P1, P2, P1nl, fnu

  fnu = 1d0 - fcb
  ! dimensionless linear power spectrum (P1)
  P1 = fc**2*CC + 2*fc*fcb*CB + fb**2*BB
  call halofit_S02(k,n,nc,knl,P1,P1nl,Om,Ov)
  P2 = fnu**2*NN + 2*fnu*fcb*(fc*CN+fb*BN)
  Pnl = P1nl*fcb**2 + P2

end subroutine halofit_nu


subroutine find_pknl_params(k,plin,rknl,rneff,rncur,nonlinear)
! * find k_NL, n_eff and n_cur (taken from CAMB)
  implicit none
  logical, intent(out) :: nonlinear
  double precision, intent(in) :: k(:), plin(:)
  double precision, intent(out) :: rknl, rneff, rncur
  double precision :: xlogr1, xlogr2, rmid, sig, d1, d2, diff

  xlogr1 = -2.0
  xlogr2 = 3.5
  nonlinear = .True.
  do
    rmid = (xlogr2+xlogr1)/2.0
    rmid = 10**rmid
    call wint(k,plin,rmid,sig,d1,d2)
    diff = sig-1.0
    if (abs(diff).le.0.001) then
      rknl  = 1./rmid
      rneff = -3-d1
      rncur = -d2   
      exit
    else if (diff.gt.0.001) then
      xlogr1 = log10(rmid)
    else if (diff.lt.-0.001) then
      xlogr2 = log10(rmid)
    end if
    if (xlogr2 < -1.9999) then !is still linear, exit
      nonlinear = .False.
      exit
    end if
  end do

end subroutine find_pknl_params


subroutine wint(k,plin,r,sig,d1,d2)
! * computing integral of Delta_L(k) with Gaussian smoothing (taken from CAMB)
  implicit none
  double precision, intent(in) :: k(:), plin(:), r
  double precision, intent(out) :: sig, d1, d2
  integer :: i, n, id
  double precision :: s(3), t, y, x, w(3), x2, rk, fac, anorm

  n = 3000
  s = 0d0
  anorm = 1/(2*pi**2)
  do i = 1, n
    t    = (i-0.5d0)/n
    y    = -1d0+1d0/t
    rk   = y
    id   = neighb(rk,k)
    d2   = plin(id)*(rk**3*anorm)
    x    = y*r
    x2   = x*x
    w(1) = exp(-x2)
    w(2) = 2*x2*w(1)
    w(3) = 4*x2*(1-x2)*w(1)
    fac  = d2/y/t/t
    s    = s + w*fac
  end do
  s   = s / dble(n)
  sig = sqrt(s(1))
  d1  = -s(2)/s(1)
  d2  = -s(2)*s(2)/s(1)/s(1) - s(3)/s(1)
      
end subroutine wint


subroutine NonLinRatios(plin,z,k,cp,pnl,ftype)
! * This implementation uses Halofit (Taken from CAMB)
!
! [inputs]
!   plin(z,k) --- linear matter power spectrum at each z and k
!   z(:)      --- redshift
!   k(:)      --- wavelength
!   cp        --- cosmological parameters
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in)  :: plin(:,:), z(:), k(:)
! (optional)
  character(*), intent(in), optional :: ftype

! [outputs]
  double precision, intent(out) :: pnl(:,:)

! [internal]
  character(128) :: fittype
  integer :: i, zi, knum, znum
  double precision :: Om_m, Om_v, a, sig, rknl, rneff, rncur, d1, d2, diff, xlogr1, xlogr2, rmid, dlin, dnl, w

  znum = size(z)
  knum = size(k)
  pnl  = plin

  fittype = 'T12'
  if (present(ftype)) fittype = ftype

  do zi = 1, znum
    ! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and
    ! curvature (rncur) of the power spectrum at the desired redshift, using method described in Smith et al (2002).
    a      = 1d0/(1d0+z(zi))
    Om_m   = omega_m(a,cp)
    Om_v   = omega_v(a,cp)
    xlogr1 = -2.0
    xlogr2 = 3.5
    do
      rmid = (xlogr2+xlogr1)/2.0
      rmid = 10**rmid
      call wint(k,plin(zi,:),rmid,sig,d1,d2)
      diff=sig-1.0
      if (abs(diff).le.0.001) then
        rknl  = 1./rmid
        rneff = -3-d1
        rncur = -d2
        exit
      elseif (diff.gt.0.001) then
        xlogr1=log10(rmid)
      elseif (diff.lt.-0.001) then
        xlogr2=log10(rmid)
      endif
      if (xlogr2 < -1.9999) then !is still linear, exit
        goto 101
      else if (xlogr2>3.4999) then !Totally crazy non-linear
        write(*,*) 'Error in halofit'
        goto 101
      end if
    end do
    ! now calculate power spectra for a logarithmic range of wavenumbers (rk)
    do i = 1, knum
      if (k(i) > 0.005d0) then
        ! dimension less linear power spectrum: dlin = k^3 * P(k) * 4*pi*V/(2*pi)^3
        dlin = plin(zi,i)*(k(i)**3/(2*pi**2))
        w = cp%w0+(1d0-a)*cp%wa
        if (fittype=='S02') call halofit_S02(k(i),rneff,rncur,rknl,dlin,dnl,Om_m,Om_v)
        if (fittype=='T12') call halofit_T12(k(i),rneff,rncur,rknl,dlin,dnl,Om_m,Om_v,cp%Om,w,cp%nu)
        if (dnl>dlin) pnl(zi,i) = (dnl/dlin)*plin(zi,i)
      end if
    end do

101 continue
  end do

end subroutine NonLinRatios


!//// Mass Function ////!

function dndM(rho,M,k,Pk)  result(f)
! mass function: dn/dM
  implicit none
  !I/O
  double precision, intent(in) :: rho, M
  double precision, intent(in) :: k(:), Pk(:)
  !internal
  double precision :: R, s2, nu, f

  R = (0.75d0*M/(pi*rho))**(1d0/3d0)
  s2 = pk2sigma(R,k,Pk)
  nu = 1.686d0**2/s2
  f = rho*nuf(nu)*(-ds2dlnM(R,k,Pk)/M)/s2

end function dndM


function ds2dlnM(R,k,Pk) result(f)
! derivative of sigma^2 w.r.t. logarithmic mass lnM
  implicit none
  !I/O
  double precision, intent(in), dimension(:) :: k, Pk
  double precision, intent(in) :: R
  !internal
  integer :: i
  double precision :: f, dk, kR

  f = 0d0
  do i = 1, size(k)-1
    dk = k(i+1)-k(i)
    kR = k(i)*R
    f = f + dk*k(i)**2*Pk(i)*tophat(kR)*dtophatdx(kR)*kR/1.5d0
  end do

end function ds2dlnM


function nuf(nu)  result(f)
! * Seth Tormen Mass Function
  implicit none
  double precision, intent(in) :: nu
  double precision :: A, p, n, f

  n = 0.707d0*nu
  A = 0.3222d0
  p = 0.3d0
  f = A*(1d0+1d0/n**p)*dsqrt(n*2d0)*dexp(-n/2d0)/dsqrt(pi)

end function nuf


function pk2sigma(R,k,Pk,n)  result(f)
! obtain sigma^2 of fluctuations using their P(k)
  implicit none
  !I/O
  integer, intent(in), optional :: n
  double precision, intent(in), dimension(:) :: k, Pk
  double precision, intent(in) :: R
  !internal
  integer :: i
  double precision :: f, dk, p

  f = 0d0
  p = 2d0
  do i = 1, size(k)-1
    dk = k(i+1)-k(i)
    if(present(n))  p = 2*dble(n) + 2d0
    f = f + dk*k(i)**p*Pk(i)*tophat(k(i)*R)**2/(2d0*pi**2)
  end do

end function pk2sigma


!//// filter ////!

function tophat(x) result(W)
! top-hat filter in Fourier space
  implicit none
  double precision, intent(in) :: x
  double precision :: W

  W = (3d0/x**3)*(dsin(x)-x*dcos(x))

end function tophat


function dtophatdx(x) result(dWdx)
! derivative of the top-hat filter
  implicit none
  double precision, intent(in) :: x
  double precision :: dWdx

  dWdx = 3d0*( (x**2-3d0)*dsin(x) + 3d0*x*dcos(x) )/x**4

end function dtophatdx


!//// Minkowski Functionals ////!

subroutine gaussian_mfs(s2,MFs,nu,n)
  implicit none 
  !I/O
  integer, intent(in) :: n
  double precision, intent(in) :: s2(0:1), nu(1:n)
  double precision, intent(out) :: MFs(n,0:2)
  !internal
  integer :: k, i
  double precision :: Ak, omega(0:2)

  omega(0:2) = [1d0,2d0,pi]
  do i = 1, n
    do k = 0, 2
      Ak = pi/(omega(2-k)*omega(k))*dsqrt(s2(1)/(2d0*s2(0)))**k/((2d0*pi)**((k+1)*0.5d0))
      MFs(i,k) = Ak*dexp(-nu(i)**2*0.5d0)*Her(k-1,nu(i))
    end do
  end do

end subroutine gaussian_mfs


end module cosmofunc

