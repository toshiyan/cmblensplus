!////////////////////////////////////////////////////!
! * Computing BB with flat sky approximation
!////////////////////////////////////////////////////!

module BB_flat
  use constants, only: pi
  use general, only: gauss_legendre_params, gl_initialize, gl_finalize
  implicit none

  private pi
  private gauss_legendre_params, gl_initialize, gl_finalize

contains


subroutine BBXY_FLAT(rL,eL,BB,XX,YY,weight,gln,gle)
! * integrating XX*YY*W^2: X=mod, Y=E or T
  implicit none
  !I/O
  character(*), intent(in), optional :: weight
  integer, intent(in) :: el(1:2), rL(1:2)
  integer, intent(in), optional :: gln
  double precision, intent(in) :: XX(:), YY(:)
  double precision, intent(in), optional :: gle
  double precision, intent(out) :: BB(:,:)
  !internal
  type(gauss_legendre_params) :: GL
  character(16) :: ftype = 'lensing'
  integer :: l, L1, L2, i, n
  double precision :: eps
  double precision, dimension(:), allocatable :: phi, intg, intc

  if (present(weight)) ftype = weight

  !* prepare GL quadrature
  n = rL(2)
  eps = 1d-14
  if (present(gln)) n = gln
  if (present(gle)) eps = gle
  call gl_initialize(GL,n,eps)

  allocate(phi(GL%n),intg(GL%n),intc(GL%n));  phi=0d0;  intg=0d0;  intc=0d0
  phi = pi + pi*GL%z

  !* calculate BB
  do l = eL(1), eL(2) ! for BB multipoles 
    L1 = rL(1)
    !* integral
    do n = 1, rL(2)-rL(1)
      L1 = L1 + 1
      do i = 1, GL%n
        L2 = int(dsqrt(dble(l)**2+dble(L1)**2-2*l*L1*dcos(phi(i))))
        call KBB(rL,l,L1,phi(i),XX(L1),YY(L2),intg(i),intc(i),ftype)
      end do
      BB(1,l) = BB(1,l) + dble(L1)*sum(intg*GL%w)*pi/(2*pi)**2
      BB(2,l) = BB(2,l) + dble(L1)*sum(intc*GL%w)*pi/(2*pi)**2
    end do
  end do

  deallocate(phi,intg,intc)
  call gl_finalize(GL)

end subroutine BBXY_FLAT


subroutine KBB(rL,l,L1,phi,XX,YY,Kg,Kc,ftype)
  implicit none
  !I/O
  character(*), intent(in) :: ftype
  integer, intent(in) :: rL(1:2), l, L1
  double precision, intent(in) :: phi, XX, YY
  double precision, intent(out) :: Kg, Kc
  !internal
  double precision :: aL0, aL1, aL2, sin2m, sin2p, sin21, sin22, cos21, L1l(2), L2l(2), L12(2), cos1, cos2, sin1, sin2, vL0(2), vL1(2), vL2(2)

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
    ! L1 * L, L1 x L
    L1l  = [dot_product(vL0,vL1),vL1(1)*vL0(2)-vL1(2)*vL0(1)]
    L2l  = [dot_product(vL0,vL2),vL2(1)*vL0(2)-vL2(2)*vL0(1)]
    ! L1 * L2, L1 x L2
    L12  = [dot_product(vL1,vL2),vL1(1)*vL2(2)-vL1(2)*vL2(1)]
    ! cos(phi), sin(phi)
    cos1 = L1l(1)/(aL1*aL0)
    sin1 = - L1l(2)/(aL1*aL0)
    cos2 = L2l(1)/(aL2*aL0)
    sin2 = - L2l(2)/(aL2*aL0)
    ! sin2(phi1), sin2(phi2), cos2(phi1)
    sin21 = 2d0*sin1*cos1
    sin22 = 2d0*sin2*cos2
    cos21 = 2d0*cos1**2 - 1d0
    ! sin2(2phi1-phi2)
    sin2m = 2d0*(sin21*cos2-sin2*cos21)*(cos21*cos2+sin21*sin2)
    select case (ftype)
    case ('lensing')
      Kg = (L12(1)*sin22)**2*XX*YY
      Kc = (L12(2)*sin22)**2*XX*YY
    case ('rotation')
      Kg = 4d0*(1d0-sin22**2)*XX*YY
    case ('patchytau')
      Kg = sin22**2*XX*YY
    case ('t2p0')
      Kg = sin21**2*XX*YY
      Kc = (1d0-sin21**2)*XX*YY
    case ('spinflip')
      Kg = sin2m**2*XX*YY
      Kc = (1d0-sin2m**2)*XX*YY
    end select
  end if

end subroutine KBB


end module BB_flat

