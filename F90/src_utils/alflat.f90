!////////////////////////////////////////////////////!
! * Computing Al with flat sky approximation
!////////////////////////////////////////////////////!

module norm_flat
  use constants, only: pi
  use general, only: gauss_legendre_params, gl_initialize, gl_finalize
  implicit none

  private pi
  private gauss_legendre_params, gl_initialize, gl_finalize

contains


subroutine ALXY_FLAT(est,eL,rL,Alg,Alc,fC,W1,W2,weight,gln,gle,lxcut)
! * integrating fXY*gXY
  implicit none
  !I/O
  character(*), intent(in) :: est
  character(*), intent(in), optional :: weight
  integer, intent(in) :: el(2), rL(2)
  integer, intent(in), optional :: gln, lxcut
  double precision, intent(in) :: fC(:), W1(:), W2(:)
  double precision, intent(in), optional :: gle
  double precision, intent(out) :: Alg(:), Alc(:)
  !internal
  type(gauss_legendre_params) :: GL
  character(16) :: ftype = 'lensing'
  integer :: l, L1, L2, i, n
  double precision :: eps, lx1, lx2
  double precision, dimension(:), allocatable :: phi, intg, intc

  write(*,*) 'alflat', est

  if (present(weight)) ftype = weight

  !* prepare GL quadrature
  n = rL(2)
  eps = 1d-14
  if (present(gln)) n = gln
  if (present(gle)) eps = gle
  call gl_initialize(GL,n,eps)

  allocate(phi(GL%n),intg(GL%n),intc(GL%n));  phi=0d0;  intg=0d0;  intc=0d0
  phi = pi + pi*GL%z

  !* calculate Al
  do l = eL(1), eL(2) ! for lensing multipoles 
    L1 = rL(1)
    !* integral
    do n = 1, rL(2)-rL(1)
      L1 = L1 + 1
      do i = 1, GL%n
        if (present(lxcut)) then
          lx1 = dble(L1)*dcos(phi(i))
          lx2 = dble(l) - lx1
          if (abs(lx1)<lxcut.or.abs(lx2)<lxcut) cycle
        end if
        L2 = int(dsqrt(dble(l)**2+dble(L1)**2-2*l*L1*dcos(phi(i))))
        select case(est)
        case ('TT')
          call KTT(rL,l,L1,phi(i),fC(L1),fC(L2),W1(L1),W2(L2),intg(i),intc(i),ftype)
        case ('EB','TB')
          call KEB(rL,l,L1,phi(i),fC(L1),W1(L1),W2(L2),intg(i),intc(i),ftype)
        end select
      end do
      Alg(l) = Alg(l) + dble(L1)*sum(intg*GL%w)*pi/(2*pi)**2
      Alc(l) = Alc(l) + dble(L1)*sum(intc*GL%w)*pi/(2*pi)**2
    end do
    if (Alg(l)/=0d0) Alg(l) = 1d0/Alg(l)
    if (Alc(l)/=0d0) Alc(l) = 1d0/Alc(l)
  end do

  deallocate(phi,intg,intc)
  call gl_finalize(GL)

end subroutine ALXY_FLAT


subroutine KTT(rL,l,L1,phi,fC1,fC2,W1,W2,Kg,Kc,ftype)
  implicit none
  !I/O
  character(*), intent(in) :: ftype
  integer, intent(in) :: rL(1:2), l, L1
  double precision, intent(in) :: phi, fC1, fC2, W1, W2
  double precision, intent(out) :: Kg, Kc
  !internal
  double precision :: al0, aL1, aL2, sin2, L1l(2), L2l(2), csp1(2), csp2(2), vL0(2), vL1(2), vL2(2)

  !* l vectors
  vl0 = [dble(l),0d0]
  vL1 = [dble(L1)*dcos(phi),dble(L1)*dsin(phi)]
  vL2 = vl0 - vL1
  al0 = dble(l)
  aL1 = dble(L1)
  aL2 = dsqrt(dot_product(vL2,vL2))

  Kg = 0d0
  Kc = 0d0
  if (rL(1)<=int(aL2).and.int(aL2)<=rL(2)) then 
    L1l = [dot_product(vl0,vL1),vL1(1)*vl0(2)-vL1(2)*vl0(1)]
    L2l = [dot_product(vl0,vL2),vL2(1)*vl0(2)-vL2(2)*vl0(1)]
    select case (ftype)
    case ('lensing')
      Kg = (fC1*L1l(1)+fC2*L2l(1))**2*W1*W2/2d0
      Kc = (fC1*L1l(2)+fC2*L2l(2))**2*W1*W2/2d0
    case ('patchytau')
      Kg = (fC1+fC2)**2*W1*W2/2d0
    end select
  end if

end subroutine KTT


subroutine KEB(rL,l,L1,phi,fC,W1,W2,Kg,Kc,ftype)
  implicit none
  !I/O
  character(*), intent(in) :: ftype
  integer, intent(in) :: rL(1:2), l, L1
  double precision, intent(in) :: phi, fC, W1, W2
  double precision, intent(out) :: Kg, Kc
  !internal
  double precision :: aL0, aL1, aL2, sin2m, sin2p, sin22, L1l(2), L2l(2), cos1, cos2, sin1, sin2, vL0(2), vL1(2), vL2(2)

  !* l vectors
  vL0 = [dble(L),0d0]
  vL1 = [dble(L1)*dcos(phi),dble(L1)*dsin(phi)]
  vL2 = vL0 - vL1
  aL0 = dble(L)
  aL1 = dble(L1)
  aL2 = dsqrt(dot_product(vL2,vL2))

  Kg = 0d0
  Kc = 0d0
  if (rL(1)<=int(aL2).and.int(aL2)<=rL(2)) then 
    ! L1 * L and L1 x L
    L1l  = [dot_product(vL0,vL1),vL1(1)*vL0(2)-vL1(2)*vL0(1)]
    L2l  = [dot_product(vL0,vL2),vL2(1)*vL0(2)-vL2(2)*vL0(1)]
    ! cos(phi), sin(phi)
    cos1 = L1l(1)/(aL1*aL0)
    sin1 = - L1l(2)/(aL1*aL0)
    cos2 = L2l(1)/(aL2*aL0)
    sin2 = - L2l(2)/(aL2*aL0)
    ! sin2(phi1-phi2)
    sin2m = 2d0*(sin1*cos2-sin2*cos1)*(cos1*cos2+sin1*sin2)
    ! sin2(phi1+phi2)
    sin2p = 2d0*(sin1*cos2+sin2*cos1)*(cos1*cos2-sin1*sin2)
    ! sin2(phi1)
    sin22 = 2d0*sin2*cos2
    select case (ftype)
    case ('lensing')
      Kg = (fC*L1l(1)*sin2m)**2*W1*W2
      Kc = (fC*L1l(2)*sin2m)**2*W1*W2
    case ('rotation')
      Kg = (2d0*fC)**2*(1d0-sin2m**2)*W1*W2
    case ('patchytau')
      Kg = (fC*sin2m)**2*W1*W2
    case ('t2p0')
      Kg = (fC*sin22)**2*W1*W2
      Kc = (fC)**2*(1d0-sin22**2)*W1*W2
    case ('spinflip')
      Kg = (fC*sin2p)**2*W1*W2
      Kc = (fC)**2*(1d0-sin2p**2)*W1*W2
    end select
  end if

end subroutine KEB


subroutine TBEB_FLAT(rL,eL,Alg,Alc,fTE,fEE,W1,W2,gln,gle,lxcut)
!* computing unnormalized cross-power spectrum TBxEB
  implicit none
  !I/O
  integer, intent(in) :: el(1:2), rL(1:2)
  integer, intent(in), optional :: gln, lxcut
  double precision, intent(in) :: fTE(:), fEE(:), W1(:), W2(:)
  double precision, intent(in), optional :: gle
  double precision, intent(out) :: Alg(:), Alc(:)
  !internal
  type(gauss_legendre_params) :: GL
  integer :: l, L1, L2, i, n
  double precision :: eps, lx1, lx2
  double precision, dimension(:), allocatable :: phi, intg, intc

  write(*,*) 'TBEB flat'

  !* prepare GL quadrature
  n = rL(2)
  eps = 1d-14
  if (present(gln)) n = gln
  if (present(gle)) eps = gle
  call gl_initialize(GL,n,eps)

  allocate(phi(GL%n),intg(GL%n),intc(GL%n));  phi=0d0;  intg=0d0;  intc=0d0
  phi = pi + pi*GL%z

  !* calculate Al
  do l = eL(1), eL(2) ! for lensing multipoles 
    L1 = rL(1)
    !* integral
    do n = 1, rL(2)-rL(1)
      L1 = L1 + 1
      do i = 1, GL%n
        if (present(lxcut)) then
          lx1 = dble(L1)*dcos(phi(i))
          lx2 = dble(l) - lx1
          if (abs(lx1)<lxcut.or.abs(lx2)<lxcut) cycle
        end if
        L2 = int(dsqrt(dble(l)**2+dble(L1)**2-2*l*L1*dcos(phi(i))))
        call KTBEB(rL,l,L1,phi(i),fTE(L1),fEE(L1),[W1(L1),W2(L2)],intg(i),intc(i))
      end do
      Alg(l) = Alg(l) + dble(L1)*sum(intg*GL%w)*pi/(2*pi)**2
      Alc(l) = Alc(l) + dble(L1)*sum(intc*GL%w)*pi/(2*pi)**2
    end do
  end do

  deallocate(phi,intg,intc)
  call gl_finalize(GL)

end subroutine TBEB_FLAT


subroutine KTBEB(rL,l,L1,phi,fTE,fEE,W,Kg,Kc)
  implicit none
  !I/O
  integer, intent(in) :: rL(1:2), l, L1
  double precision, intent(in) :: phi, fTE, fEE, W(2)
  double precision, intent(out) :: Kg, Kc
  !internal
  double precision :: al0, aL1, aL2, sin2, L1l(2), L2l(2), csp1(2), csp2(2), vL0(2), vL1(2), vL2(2)

  !* l vectors
  vl0 = [dble(l),0d0]
  vL1 = [dble(L1)*dcos(phi),dble(L1)*dsin(phi)]
  vL2 = vl0 - vL1
  al0 = dble(l)
  aL1 = dble(L1)
  aL2 = dsqrt(dot_product(vL2,vL2))

  Kg = 0d0
  Kc = 0d0
  if (rL(1)<=int(aL2).and.int(aL2)<=rL(2)) then 
    L1l = [dot_product(vl0,vL1),vL1(1)*vl0(2)-vL1(2)*vl0(1)]
    L2l = [dot_product(vl0,vL2),vL2(1)*vl0(2)-vL2(2)*vl0(1)]
    csp1 = [L1l(1)/(aL1*al0),- L1l(2)/(aL1*al0)]
    csp2 = [L2l(1)/(aL2*al0),- L2l(2)/(aL2*al0)]
    sin2 = 2d0*(csp1(2)*csp2(1)-csp2(2)*csp1(1))*(csp1(1)*csp2(1)+csp1(2)*csp2(2))
    Kg = (fTE*L1l(1)*sin2)*(fEE*L1l(1)*sin2)*W(1)*W(2)
    Kc = (fTE*L1l(2)*sin2)*(fEE*L1l(2)*sin2)*W(1)*W(2)
  end if

end subroutine KTBEB


end module norm_flat


#ifdef old

!////////////////////////////////////////////////////!
! * Old nldd_flat code
!////////////////////////////////////////////////////!

subroutine AL_FLAT(GL,O,Q,E,el,tl,ucl,lcl,Al,Nl)
  implicit none
  !I/O
  type(GAUSS_LEGENDRE_PARAMS), intent(in) :: GL
  type(OBSERVABLES), intent(in) :: O
  type(QUADRATIC_COMBINATION), intent(in) :: Q
  type(LENSING_ESTIMATOR), intent(in) :: E
  integer, intent(in) :: el(2), tl(2)
  real(dl), intent(in) :: ucl(:,:), lcl(:,:)
  real(dl), intent(inout) :: Al(:,:,:)
  real(dl), intent(out) :: Nl(:,:,:)
  !internal
  integer :: jmax, j, l, k, i, n, a, b
  real(dl) :: LL, dk, phi
  real(dl), dimension(:,:,:), allocatable :: Temp, dAL, Ij
  real(dl), dimension(:,:,:,:), allocatable :: Aj

  jmax = el(2)-el(1)+1

  !#### calculate noise power spectrum ####
  allocate(Aj(Q%MV,E%n,E%n,jmax))
  Aj = 0d0
  Nl = 0d0
  !For Gauss-Legendre Integration
  allocate(Temp(Q%MV,E%n,E%n),dAl(Q%MV,E%n,E%n))
  do j = 1, jmax ! for lensing multipoles 
    l = el(1) + j - 1
    write(*,*) "l=", l
    k = tl(1)
    dk = 1
    Aj(:,:,:,j) = 0d0

    !* integral *!
    !the lensing multipole vector is fixed to x-axis since the results 
    !are independent of the direction of lensing multipole
    allocate(Ij(4,E%n,E%n))
    do n = 1, tL(2)-tL(1)
      k = k + dk
      Temp = 0d0 !note: xx=0, xm=pi, xl=pi
      do i = 1, GL%n
        phi = pi + pi*GL%z(i)
        LL = dsqrt(l**2 + k**2 - 2*l*k*dcos(phi))
        if (tL(1)<=LL.and.LL<=tL(2)) then 
          call al_kernel_flat(O,Q,E,l,k,phi,el(2),ucl,lcl,dAl,Ij)
          Temp = Temp + GL%w(i)*pi*dAl
        end if
      end do
      Aj(:,:,:,j) = Aj(:,:,:,j) + dble(k*dk)*Temp/(2*pi)**2
    end do

    do i = 1, Q%n
      n = 1
      do a = 1, E%n
        !Normalization (a = b)
        if(Aj(i,a,a,j)/=0d0) Al(i,n,l) = 1d0/Aj(i,a,a,j)
        !Response function (a not b)
        do b = a + 1, E%n
          if(Aj(i,a,a,j)/=0d0) Al(i,n,l) = Aj(i,a,b,j)/Aj(i,a,a,j)
          n = n + 1
        end do
      end do
    end do
    !need fix
    !if(Q%n/=1) call OPTCOMB(Q,Al(:,:,l),Ij(:,:,1),Nl(:,:,l)) !MV
    deallocate(Ij)
  end do

  deallocate(Temp,dAl)

end subroutine AL_FLAT


subroutine AL_KERNEL_FLAT(O,Q,E,ell,L1,tp,iLmax,FC,OC,K,Il)
  implicit none
  !I/O
  type(OBSERVABLES), intent(in) :: O
  type(QUADRATIC_COMBINATION), intent(in) :: Q
  type(LENSING_ESTIMATOR), intent(in) :: E
  integer, intent(in) :: ell, L1, iLmax
  real(dl), intent(in) :: tp
  real(dl), intent(in) :: FC(:,:), OC(:,:)
  real(dl), intent(out) :: K(:,:,:), Il(:,:,:)
  !internal
  integer :: L2, XY, a, b
  integer, parameter :: x=1, y=2
  real(dl) :: L1gl,L2gl,L1cl,L2cl 
  real(dl) :: cospL1, sinpL1, cospL2, sinpL2
  real(dl) :: cos2pL1, sin2pL1, cos2pL2, sin2pL2
  real(dl) :: cos2pL1L2, sin2pL1L2
  real(dl) :: al0, aL1, aL2
  real(dl) :: vl0(2), vL1(2), vL2(2)
  real(dl) :: A_l1l2, A_l2l1, B_l1l2, C_l1l2
  real(dl) :: f(E%n,Q%n)
  real(dl) :: ggTE_l1l2, ggTE_l2l1, fgTE_l1l2, fgTE_l2l1
  real(dl) :: gcTE_l1l2, gcTE_l2l1, fcTE_l1l2, fcTE_l2l1
  real(dl) :: gmTE_l1l2, gmTE_l2l1, fmTE_l1l2, fmTE_l2l1

  f = 0d0
  al0 = dble(ell)
  aL1 = dble(L1)

  !* l-space
  vl0(x) = al0
  vl0(y) = 0d0
  vL1(x) = aL1*dcos(tp)
  vL1(y) = aL1*dsin(tp)
  vL2 = vl0 - vL1

  !* inner dot products 
  L1gl = dot_product(vl0,vL1)
  L2gl = dot_product(vl0,vL2)
  L1cl = vL1(x)*vl0(y) - vL1(y)*vl0(x)
  L2cl = vL2(x)*vl0(y) - vL2(y)*vl0(x)

  !* sin and cosin
  aL2 = dsqrt(dot_product(vL2,vL2))
  L2 = int(aL2)
  cospL1 = L1gl/(aL1*al0)
  sinpL1 = - L1cl/(aL1*al0)
  cospL2 = L2gl/(aL2*al0)
  sinpL2 = - L2cl/(aL2*al0)
  cos2pL1 = 2*cospL1**2 - 1.d0
  sin2pL1 = 2*cospL1*sinpL1
  cos2pL2 = 2*cospL2**2 - 1.d0
  sin2pL2 = 2*cospL2*sinpL2
  cos2pL1L2 = cos2pL1*cos2pL2 + sin2pL1*sin2pL2
  sin2pL1L2 = sin2pL1*cos2pL2 - cos2pL1*sin2pL2

  !* weight functions
  if(E%DO(E%g)) then
    f(E%g,Q%TT) = FC(O%TT,L1)*L1gl + FC(O%TT,L2)*L2gl
    f(E%g,Q%TB) = FC(O%TE,L1)*L1gl*sin2pL1L2
    f(E%g,Q%EE) = (FC(O%EE,L1)*L1gl + FC(O%EE,L2)*L2gl)*cos2pL1L2
    f(E%g,Q%EB) = (FC(O%EE,L1)*L1gl + FC(O%BB,L2)*L2gl)*sin2pL1L2
    f(E%g,Q%BB) = (FC(O%BB,L1)*L1gl + FC(O%BB,L2)*L2gl)*cos2pL1L2
  end if
  if(E%DO(E%c)) then
    f(E%c,Q%TT) = FC(O%TT,L1)*L1cl + FC(O%TT,L2)*L2cl
    f(E%c,Q%TB) = FC(O%TE,L1)*L1cl*sin2pL1L2
    f(E%c,Q%EE) = (FC(O%EE,L1)*L1cl + FC(O%EE,L2)*L2cl)*cos2pL1L2
    f(E%c,Q%EB) = (FC(O%EE,L1)*L1cl + FC(O%BB,L2)*L2cl)*sin2pL1L2
    f(E%c,Q%BB) = (FC(O%BB,L1)*L1cl + FC(O%BB,L2)*L2cl)*cos2pL1L2
  end if
  if(E%DO(E%m)) then
    f(E%m,Q%TT) = FC(O%TT,L1) + FC(O%TT,L2)
    !if(pEB) then
    !  f(E%m,QTB) = 0.d0
    !  f(E%m,QEE) = (aL1/aL2)**2*FC(O%EE,L1) + (aL2/aL1)**2*FC(O%EE,L2)
    !  f(E%m,QEB) = 0.d0
    !  f(E%m,QBB) = (aL1/aL2)**2*FC(O%BB,L1) + (aL2/aL1)**2*FC(O%BB,L2)
    f(E%m,Q%TB) = FC(O%TE,L1)*sin2pL1L2
    f(E%m,Q%EE) = (FC(O%EE,L1) + FC(O%EE,L2))*cos2pL1L2
    f(E%m,Q%EB) = (FC(O%EE,L1) + FC(O%BB,L2))*sin2pL1L2
    f(E%m,Q%BB) = (FC(O%BB,L1) + FC(O%BB,L2))*cos2pL1L2
  end if
  if(E%DO(E%t)) then
    f(E%t,Q%EB) = -2*FC(O%TE,L1)*sin2pL2
  end if
  if(E%DO(E%p)) then
    f(E%p,Q%EB) = 2*(FC(O%EE,L2)-FC(O%BB,L1))*(sin2pL1*cos2pL2+sin2pL2*cos2pL1)
  end if
  !TE
  fgTE_l1l2 = FC(O%TE,L1)*L1gl*cos2pL1L2 + FC(O%TE,L2)*L2gl
  fgTE_l2l1 = FC(O%TE,L2)*L2gl*cos2pL1L2 + FC(O%TE,L1)*L1cl
  fcTE_l1l2 = FC(O%TE,L1)*L1cl*cos2pL1L2 + FC(O%TE,L2)*L2cl
  fcTE_l2l1 = FC(O%TE,L2)*L2cl*cos2pL1L2 + FC(O%TE,L1)*L1cl
  fmTE_l1l2 = FC(O%TE,L1)*cos2pL1L2 + FC(O%TE,L2)
  fmTE_l2l1 = FC(O%TE,L2)*cos2pL1L2 + FC(O%TE,L1)
  A_l1l2 = OC(O%TT,L2)*OC(O%EE,L1)
  A_l2l1 = OC(O%TT,L1)*OC(O%EE,L2)
  B_l1l2 = OC(O%TE,L1)*OC(O%TE,L2)
  C_l1l2 = OC(O%TT,L1)*OC(O%EE,L2)*OC(O%TT,L2)*OC(O%EE,L1) - OC(O%TE,L1)**2*OC(O%TE,L2)**2
  ggTE_l1l2 = (A_l1l2*fgTE_l1l2 - B_l1l2*fgTE_l2l1)/C_l1l2
  gcTE_l1l2 = (A_l1l2*fcTE_l1l2 - B_l1l2*fcTE_l2l1)/C_l1l2
  gmTE_l1l2 = (A_l1l2*fmTE_l1l2 - B_l1l2*fmTE_l2l1)/C_l1l2

  !**** need check ****!
  Il = 0d0
  !do a = 1, emax
  !  Il(Q%TTTE,a,a) = f(a,TT)/(2.d0*OC(TT,lp)*OC(TT,lpp))*(ggTE_l1l2*OC(TT,lp)*OC(TE,lpp)+ggTE_l2l1*OC(TE,lp)*OC(TT,lpp))
  !  Il(Q%TTEE,a,a) = f(a,TT)*2._dl*ggEE*OC(TE,lp)*OC(TE,lpp)
  !  Il(Q%TEEE,a,a) = ggTE_l1l2*2._dl*ggEE*OC(TE,lp)*OC(EE,lpp)
  !  Il(Q%TBEB,a,a) = ggTB*ggEB*OC(TE,lp)*OC(BB,lpp)*sq_sin2phi
  !end do
  !********************!

  do XY = 1, Q%n
    if(XY==Q%TE) then 
      if(E%CLDO(E%gg)) K(Q%TE,E%g,E%g) = ggTE_l1l2*fgTE_l1l2
      if(E%CLDO(E%gm)) K(Q%TE,E%g,E%m) = ggTE_l1l2*fmTE_l1l2
      if(E%CLDO(E%cc)) K(Q%TE,E%c,E%c) = gcTE_l1l2*fcTE_l1l2
      if(E%CLDO(E%mm)) K(Q%TE,E%m,E%m) = gmTE_l1l2*fmTE_l1l2
      !K(Q%TE,E%t,E%t) = 0.d0!f2TE_l1l2*g2TE_l1l2
      !K(Q%TE,E%p,E%p) = 0.d0!f4TE_l1l2*g4TE_l1l2
      !K(Q%TE,E%g,E%p) = 0.d0!f2TE_l1l2*g2TE_l1l2
      !K(Q%TE,E%g,E%p) = 0.d0!f4TE_l1l2*g4TE_l1l2
    else 
      do a = 1, E%n
        do b = a, E%n
          K(XY,a,b) = integrand(O,Q,XY,f(a,XY),f(b,XY),OC(:,L1),OC(:,L2))
        end do
      end do
    end if
  end do


end subroutine AL_KERNEL_FLAT


function integrand(O,Q,XY,fa,fb,OC1,OC2)
  implicit none 
  !I/O
  type(OBSERVABLES), intent(in) :: O
  type(QUADRATIC_COMBINATION), intent(in) :: Q
  integer, intent(in) :: XY
  real(dl), intent(in) :: fa, fb, OC1(:), OC2(:)
  real(dl) :: integrand, g
  !internal

  if(XY==Q%TT) then 
    g = fb/(2.d0*OC1(O%TT)*OC2(O%TT))
  else if(XY==Q%EE) then 
    g = fb/(2.d0*OC1(O%EE)*OC2(O%EE))
  else if(XY==Q%BB) then 
    g = fb/(2.d0*OC1(O%BB)*OC2(O%BB))
  else if(XY==Q%TB) then 
    g = fb/(OC1(O%TT)*OC2(O%BB))
  else if(XY==Q%EB) then 
    g = fb/(OC1(O%EE)*OC2(O%BB))
  end if

  integrand = fa*g

end function integrand 


#endif old

