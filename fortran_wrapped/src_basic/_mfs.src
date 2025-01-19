!////////////////////////////////////////////////
! Kernel of Analytic MFs upto 4th order
!////////////////////////////////////////////////

module mfs_kernel
  use myconst, only: dl, pi
  use myutils 
  use funcs, only: Her, w00_ini, w3j_recursion
  implicit none

  private dl, pi
  private Her, W00_ini, w3j_recursion
  
contains


subroutine CalcMFs(sigma2,Skew,Kurt,MFs,nus)
  implicit none 
  !I/O
  double precision, intent(in) :: sigma2(0:1), Skew(:), Kurt(:), nus(:)
  double precision, intent(out) :: MFs(:,:)
  !internal
  integer :: k, i
  double precision :: vk0, vk1, vk2, Ak, nu

  do i = 1, size(MFs,dim=2)
    do k = 0, 2
      nu = nus(i)
      vk0 = Her(k-1,nu)
      vk1 = Her(k+2,nu)*Skew(1)/6.d0 - Her(k,nu)*Skew(2)*dble(k)*0.25d0 - Her(k-2,nu)*Skew(3)*dble(k*(k-1))*0.25d0
      select case(k)
      case(0)
        vk2 = skew(1)**2/72.d0*Her(5,nu) + Kurt(1)/24.d0*Her(3,nu)
      case(1)
        vk2 = skew(1)**2/72.d0*Her(6,nu) + (Kurt(1)-skew(1)*skew(2))/24.d0*Her(4,nu) - (Kurt(2)+3.d0*skew(2)**2/8.d0)*Her(2,nu)/12.d0 - Kurt(4)/8.d0
      case(2)
        vk2 = skew(1)**2/72.d0*Her(7,nu) + (Kurt(1)-2*skew(1)*skew(2))/24.d0*Her(5,nu) - (Kurt(2)+skew(2)**2/2.d0)*Her(3,nu)/6.d0 - (Kurt(3)+skew(2)*skew(3)/2.d0)*Her(1,nu)/2.d0
      end select
      Ak = pi/(omega(2-k)*omega(k))*dsqrt(sigma2(1)/(2.d0*sigma2(0)))**k/((2.d0*pi)**((k+1)*0.5d0))
      MFs(k+1,i) = Ak*dexp(-nu**2*0.5d0)*( vk0 + dsqrt(sigma2(0))*vk1 + sigma2(0)*vk2 )
    end do
  end do

end subroutine CalcMFs


function omega(n)
  implicit none
  integer, intent(in) :: n 
  double precision :: omega

  select case (n)
  case (0)
    omega = 1d0
  case (1)
    omega = 2d0
  case default
    omega = pi
  end select

end function omega


subroutine calc_kurtosis(K,sigma2,CTT,Cdd,eL,th)
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: sigma2(0:1),CTT(:),Cdd(:),th
  double precision, intent(inout) :: K(:)
  !internal
  double precision :: W_sup, W_mid, W_inf 
  integer :: L, L1, L2, L2max, n, nmax
  double precision :: Lc, q2, mcA(3), fac

  K = 0d0
  q2 = 0.5d0*sigma2(1)/sigma2(0)

  do L=2, eL(2)
    mcA = 0d0
    write(*,*) L
    Lc = eL(1) + 1
    if (dble(L)*0.5d0 >= eL(1)+1) Lc = dble(L)*0.5d0
    if (dint(Lc).ne.Lc) Lc = Lc+0.5d0
    if (Lc.le.2d0) Lc = 2d0

    !* In what follows, loop is run over the case; L1+1 <= L2 
    !    (i.e., L1 >= L2 is added to L1+1 <= L2). 
    !* The range of L2 is restricted by the triangle relation: 
    !       |L-L1| <= L2 <= L+L1 
    !* L2min = max(L1,L-L1,tlmin)
    !    If L1 < L/2, L2min = max(L-L1,tlmin) > L1
    !    If L1 >= L/2, L2min = L1+1 >= L-L1

    do L1 = eL(1), eL(2)
      L2max = L + L1
      !* case L2 = L2max
      L2 = L2max
      W_sup = w00_ini(L,L1)
      if (raneq(L2,eL)) call calc_mcA(mcA,L,L1,L2,q2,W_sup,CTT,eL,th)
      !* case L2 = L2max - 1
      L2 = L2max - 1
      W_mid = 0d0
      if (raneq(L2,eL)) call calc_mcA(mcA,L,L1,L2,q2,W_mid,CTT,eL,th)
      !* sum over other cases, down to max (L2 = L - L1, tlmin)
      if (L1<int(Lc)) nmax = min(2*L1-1, L+L1-1-eL(1)) 
      !* sum over other cases, down to L2 = L1+1 ( > L/2 )
      if (L1>=int(Lc)) nmax = L-1
      do n = 1, nmax
        L2 = L2max-1-n
        call Recursion(L,L1,L2,0,0,W_sup,W_mid,W_inf)
        if (raneq(L2,eL)) call calc_mcA(mcA,L,L1,L2,q2,W_inf,CTT,eL,th)
        W_sup  = W_mid
        W_mid  = W_inf
      end do
    end do

    fac = Cdd(L)/dble(2*L+1)
    K(1) = K(1) + fac*mcA(1)**2
    K(2) = K(2) - fac*mcA(1)*mcA(2)*0.25d0/q2
    K(3) = K(3) + fac*(mcA(3)**2-mcA(2)**2)*(0.25d0/q2)**2
    K(4) = K(4) + fac*(mcA(2)-mcA(3))**2*0.5d0*(0.25d0/q2)**2
  end do

  K = K*0.25d0/(pi*sigma2(0)**3)
  write(*,*) "kurtosis", K(1:4)

end subroutine calc_kurtosis


subroutine calc_mcA(mcA,L,L1,L2,q2,W,CLTT,eL,th)
  implicit none
  !I/O
  double precision, intent(inout) :: mcA(3)
  integer, intent(in) :: L, L1, L2, eL(2)
  double precision, intent(in) :: q2, W, CLTT(:), th
  !internal
  double precision :: F, I, fac_12, fac_21
  double precision :: aL, aL1, aL2

  aL = dble(L)
  aL1 = dble(L1)
  aL2 = dble(L2)
  fac_12 = aL*(aL+1) + aL1*(aL1+1) - aL2*(aL2+1)
  fac_21 = aL*(aL+1) + aL2*(aL2+1) - aL1*(aL1+1)
  I = W*dsqrt((2*aL+1)*(2*aL1+1)*(2*aL2+1)/(4.d0*pi))
  F = I**2*Wl(eL,L1,th)*Wl(eL,L2,th)*(CLTT(L1)*fac_12+CLTT(L2)*fac_21)
  if(L1==L2) F = F*0.5d0

  mcA(1) = mcA(1) + F
  mcA(2) = mcA(2) + F*( aL1*(aL1+1) + aL2*(aL2+1) )
  mcA(3) = mcA(3) + F*aL*(aL+1)

end subroutine calc_mcA


subroutine Recursion(L1,ell,L2,M1,M2,W_sup,W_mid,W_inf)
  implicit none
  integer, intent(in) :: ell, L1, L2, M1, M2
  double precision, intent(in) :: W_sup, W_mid
  double precision, intent(out) :: W_inf
  double precision :: W3j(3)

  W3j(1) = W_sup
  W3j(2) = W_mid
  call w3j_recursion(L1,ell,L2,M1,M2,W3j)
  W_inf = W3j(3)

end subroutine Recursion


subroutine calc_skew(Skew,sigma2,CTT,CTd,eL,th)
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(inout) :: Skew(:)
  double precision, intent(in) :: sigma2(0:1), CTT(:), CTd(:), th
  !internal
  logical :: even
  double precision :: W_sup, W_mid, W_inf 
  integer L1, L2, L3, L2max, nn 
  double precision :: aL1, aL2, aL3, mcS(3), q2
  double precision :: a, c, iniW, iniW0, bW, WL1L2p0, WL1L2p1

  q2 = 0.5_dl*sigma2(1)/sigma2(0)

  ! * range of L0 --- [2, lmax] 
  ! * range of L1 --- [L0, lmax] 
  ! * range of L2 --- [L1, lmax] 
  !     The range of L2 is restricted by the triangle relation:
  !         |L0-L1| <= L2 <= L0+L1 
  ! * finally multiplied by "6" 
  ! * initial wigner-3j symbols are recursively computed: 
  !   (2,2,2) - (3,3,4) - (4,4,4) - (5,5,6) - ...

  iniW0 = 1.d0
  do L1 = 2, eL(2)
    aL1 = dble(L1)
    call Init_Wjm(L1,iniW,iniW0,even)
    if(even) then !for L0 = 2, 4, 6, ...
      WL1L2p0 = iniW
      WL1L2p1 = 0.d0
    else !for L0 = 3, 5, 7, ...
      WL1L2p0 = 0.d0
      WL1L2p1 = iniW
    end if
    do L2 = L1, eL(2)
      aL2 = dble(L2)
      if(even) then !L1, L2, L2
        L3 = L2
        aL3 = aL2
        call reducedbisp(L1,L2,L2,bW,CTT,CTd,eL,th)
        call calc_prefactor(mcS,aL1,aL2,aL3,q2)
        W_inf = WL1L2p0
        W_mid = 0.d0
        skew = skew + mcS*bW*W_inf**2
      else
        if(L2+1<=eL(2)) then
          L3 = L2 + 1
          aL3 = aL2+1
          call reducedbisp(L1,L2,L2+1,bW,CTT,CTd,eL,th)
          call calc_prefactor(mcS,aL1,aL2,aL3,q2)
          W_inf = WL1L2p0
          W_mid = WL1L2p1
          skew = skew + mcS*bW*W_inf**2
        end if 
      end if 
      if (L2+2<=eL(2)) then
        do L3 = L2+2, min(eL(2),L1+L2)
          aL3 = dble(L3)
          ! recursion formula 
          call Recursion_Ini(aL1,aL2,aL3,a,c)
          call reducedbisp(L1,L2,L3,bW,CTT,CTd,eL,th)
          call calc_prefactor(mcS,aL1,aL2,aL3,q2)
          W_sup = -W_inf*(aL3-1.d0)*a/(aL3*c) 
          skew = skew + mcS*bW*W_inf**2
          ! shift for the next recursion
          W_inf  = W_mid
          W_mid  = W_sup
        end do 
      end if
      WL1L2p0 = -dsqrt((2*aL2-aL1+2)*(2*aL2-aL1+1)/(2*aL2+aL1+3)/(2*aL2+aL1+2))*((2*aL2+aL1+2)/(2*aL2-aL1+2))*WL1L2p0
      WL1L2p1 = -dsqrt((2*aL2-aL1+3)*(2*aL2-aL1+2)/(2*aL2+aL1+4)/(2*aL2+aL1+3))*((2*aL2+aL1+3)/(2*aL2-aL1+3))*WL1L2p1
    end do 
    iniW0 = iniW
  end do

  skew(:) = 1.5_dl/(pi*sigma2(0)**2)*skew(:)
  !skew(2) = 0.375_dl/(pi*sigma2(0)*sigma2(1))*skew(2)
  !skew(3) = 0.75_dl/(pi*sigma2(1)**2)*skew(3)
  write(*,*) "skew", skew(1), skew(2), skew(3)

end subroutine calc_skew


subroutine Init_Wjm(L0,W,W0,even)
  implicit none
  !I/O
  logical, intent(out) :: even
  integer, intent(in) :: L0
  double precision, intent(out) :: W
  double precision, intent(in) :: W0
  !internal
  integer :: a, i
  double precision :: d

  even = .false.
  a = (L0-1)/2
  if(mod(L0,2)==0) then 
    even = .true.
    a = L0/2
  end if
  d = dble(a)

  ! "even" to "odd" or "odd" to "even"
  if(even) then 
    W = -W0*dsqrt((6.d0*d-3.d0)/(6.d0*d+1.d0))
  else
    W = W0*dsqrt(6.d0*(6.d0*d+4.d0)*(3.d0*d+1.d0)*(2.d0*d+1.d0)) &
      /(d+1.d0)/dsqrt((6.d0*d+5.d0))
  end if

  !for (2,2,2)
  if(L0==2) W = - dsqrt(2.d0/35.d0)

end subroutine Init_Wjm


subroutine Recursion_Ini(L,L1,L2,a,c)
  implicit none
  !I/O
  double precision, intent(in) :: L, L1, L2
  double precision, intent(out) :: a, c

  a = L2*dsqrt((-L2+L+L1+2)*(L2-L+L1-1)*(L2+L-L1-1)*(L2+L+L1))
  c = (L2-1)*dsqrt((-L2+L+L1+1)*(L2-L+L1)*(L2+L-L1)*(L2+L+L1+1))

end subroutine Recursion_Ini


subroutine calc_prefactor(mcS,aL1,aL2,aL3,q2)
  implicit none
  double precision, intent(inout) :: mcS(3)
  double precision, intent(in) :: aL1, aL2, aL3, q2
  
  mcS(1) = 1.d0
  mcS(2) = -(aL1*(aL1+1)+aL2*(aL2+1)+aL3*(aL3+1))/6.d0/q2
  mcS(3) = ((aL1*(aL1+1)+aL2*(aL2+1)-aL3*(aL3+1))*aL3*(aL3+1) & 
    +(aL2*(aL2+1)+aL3*(aL3+1)-aL1*(aL1+1))*aL1*(aL1+1) & 
    +(aL3*(aL3+1)+aL1*(aL1+1)-aL2*(aL2+1))*aL2*(aL2+1))/12.d0/q2

end subroutine calc_prefactor


subroutine reducedbisp(L1,L2,L3,bW,CTT,CTd,eL,th)
  implicit none 
  !I/O
  integer, intent(in) :: L1, L2, L3, eL(2)
  double precision, intent(in) :: CTT(:), CTd(:), th
  double precision, intent(out) :: bW
  !internal
  double precision :: b123, b132, b213, b231, b312, b321, aL1, aL2, aL3
  double precision :: a2L1, a2L2, a2L3   

  aL1 = dble(L1)
  aL2 = dble(L2)
  aL3 = dble(L3)
  a2L1 = aL1*(aL1+1.d0)
  a2L2 = aL2*(aL2+1.d0)
  a2L3 = aL3*(aL3+1.d0)
  b123 = 0.5d0*(a2L2+a2L3-a2L1)*CTd(L2)*CTT(L3)
  b132 = 0.5d0*(a2L3+a2L2-a2L1)*CTd(L3)*CTT(L2)
  b213 = 0.5d0*(a2L1+a2L3-a2L2)*CTd(L1)*CTT(L3)
  b231 = 0.5d0*(a2L3+a2L1-a2L2)*CTd(L3)*CTT(L1)
  b312 = 0.5d0*(a2L1+a2L2-a2L3)*CTd(L1)*CTT(L2)
  b321 = 0.5d0*(a2L2+a2L1-a2L3)*CTd(L2)*CTT(L1)
  bW = b123 + b132 + b213 + b231 + b312 + b321
  bW = bW*Wl(eL,L1,th)*Wl(eL,L2,th)*Wl(eL,L3,th) &
    *(2.d0*aL1+1.d0)*(2.d0*aL2+1.d0)*(2.d0*aL3+1.d0)/(4.d0*pi)

end subroutine reducedbisp


function WL(eL,L,th)
  implicit none
  integer, intent(in) :: L, eL(2)
  double precision :: th
  double precision :: WL

  if(L>=eL(1).and.L<=eL(2)) then
    WL = dexp(-dble(L*(L+1))*th**2*0.5d0)
  else 
    WL = 0.d0
  end if

end function WL


function WLL(eL,th,lmax)
  implicit none
  integer, intent(in) :: lmax, eL(2)
  integer :: l
  double precision :: th
  double precision :: WLL(lmax)

  do L = 1, lmax
    WLL(L) = WL(eL,L,th)
  end do

end function WLL


end module mfs_kernel

