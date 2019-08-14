!///////////////////////////////////////////////////////////////////////!
! Lensing Bispectrum
!///////////////////////////////////////////////////////////////////////!

module bispec
  use general,   only: GLdxs, GLpoints, linspace, gauss_legendre_params, gl_initialize, gl_finalize
  use cosmofunc, only: cosmoparams, set_cosmoparams, omega_m
  use pstool,    only: binned_ells
  use bstool,    only: bispecfunc, bispec_lens_bin, bispec_lens_lss, bispec_lens_pb, bispec_lens_lss_init, bispec_lens_pb_init, zinterp, skewspec_lens
  implicit none

contains


subroutine bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,bl0,bl1,pktype,ltype)
  implicit none
  !I/O
  character(8), intent(in) :: shap, cpmodel, model
  integer, intent(in) :: lmin, lmax, zn, kn
  double precision, intent(in) :: lan, kan, zs
  double precision, intent(in), dimension(1:zn) :: z, dz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(out), dimension(lmin:lmax) :: bl0, bl1
  !optional
  character(4), intent(in), optional :: pktype, ltype
  !f2py character(4) :: pktype = 'T12'
  !f2py character(4) :: ltype  = ''
  !internal
  character(4) :: pkfit = 'T12', lt = ''
  integer :: l0, oL(2), i
  double precision, allocatable :: wp(:,:), ck(:,:)
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  ! cosmological parameters (to compute background quantities)
  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  if (lmin<1) stop 'lmin should be >=1' 
  oL = (/lmin,lmax/)

  ! precompute quantities for bispectrum
  allocate(wp(zn,lmax),ck(zn,lmax))
  if (present(pktype)) pkfit = pktype
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,oL,model,pkfit) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,oL,wp,ck)

  ! MG parameter z-evolution
  allocate(b%mgp(2,zn));  b%mgp=1d0
  do i = 1, zn
    if (lan/=0d0) b%mgp(1,i) = (omega_m(1d0/(1d0+z(i)),cp))**lan
    if (kan/=0d0) b%mgp(2,i) = (omega_m(1d0/(1d0+z(i)),cp))**kan
  end do
  
  ! bispectrum (density and post-Bron bispectrum)
  write(*,*) 'compute bispectrum', shap
  l0 = lmin+(lmax-lmin)*0.5d0/20d0
  if (mod(l0,2)/=0) l0 = l0+1
  if (present(ltype)) lt = ltype
  call bispec_lens_lss(cp,b,shap,oL,model,ltype,l0,bl0)
  call bispec_lens_pb(shap,oL,wp,ck,bl1,l0,ltype)

end subroutine bispeclens


subroutine bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,bc,bl0,bl1,pktype,ltype)
  implicit none
  !I/O
  character(8), intent(in) :: shap, cpmodel, model
  integer, intent(in) :: lmin, lmax, zn, bn, kn
  double precision, intent(in) :: lan, kan, zs
  double precision, intent(in),  dimension(1:zn) :: z, dz
  double precision, intent(in),  dimension(1:kn) :: k, pk0
  double precision, intent(out), dimension(1:bn) :: bc, bl0, bl1
  !optional
  character(4), intent(in), optional :: pktype, ltype
  !f2py character(4) :: pktype = 'T12'
  !f2py character(4) :: ltype  = ''
  !internal
  character(4) :: pkfit = 'T12', lt = ''
  integer :: oL(2), l, i, eL1(2), eL2(2), eL3(2)
  double precision, allocatable :: bp(:), wp(:,:), ck(:,:)
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  if (lmin<1) stop 'lmin should be >=1' 
  oL = (/lmin,lmax/)

  ! precompute quantities for bispectrum
  allocate(wp(zn,oL(2)),ck(zn,oL(2)))
  if (present(pktype)) pkfit = pktype
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,oL,model,pkfit) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,oL,wp,ck)

  ! MG parameter z-evolution
  allocate(b%mgp(2,zn));  b%mgp=1d0
  do i = 1, zn
    if (lan/=0d0) b%mgp(1,i) = (omega_m(1d0/(1d0+z(i)),cp))**lan
    if (kan/=0d0) b%mgp(2,i) = (omega_m(1d0/(1d0+z(i)),cp))**kan
  end do
  
  ! binning
  allocate(bp(bn+1))
  call binned_ells(oL,bp,bc)

  ! bispectrum (density and post-Bron bispectrum)
  write(*,*) 'compute bispectrum', shap
  if (present(ltype)) lt = ltype
  do l = 1, bn
    select case(shap)
    case ('equi')
      eL1 = int(bp(l:l+1))
      eL2 = int(bp(l:l+1))
      eL3 = int(bp(l:l+1))
    case ('fold')
      eL1 = int(bp(l:l+1))
      eL2 = [max(2,int(bp(l)/2)),int(bp(l+1)/2)]
      eL3 = [max(2,int(bp(l)/2)),int(bp(l+1)/2)]
    case ('sque')
      eL1 = int(bp(1:2))
      eL2 = int(bp(l:l+1))
      eL3 = int(bp(l:l+1))
    case ('angl')
      eL1 = int(bp(l:l+1))
      eL2 = int(bp(bn/2:bn/2+1))
      eL3 = int(bp(bn/2:bn/2+1))
    end select
    write(*,*) eL1, eL2, eL3
    call bispec_lens_bin(cp,b,eL1,eL2,eL3,wp,ck,model,bl0(l),bl1(l))
  end do

end subroutine bispeclens_bin


subroutine zpoints(zmin,zmax,zn,z,dz,zspace)
  implicit none
  !I/O
  integer, intent(in) :: zn
  double precision, intent(in) :: zmin, zmax
  double precision, intent(out), dimension(1:zn) :: z, dz
  integer, intent(in), optional :: zspace
  !f2py integer :: zspace = 1

  ! precomputing interpolation points for z
  if (present(zspace)) then
    call zinterp(zmin,zmax,zn,zspace,z,dz)
  else
    call zinterp(zmin,zmax,zn,1,z,dz)
  end if

end subroutine zpoints


subroutine skewspeclens(cpmodel,model,z,dz,zn,zs,olmin,olmax,lmin,lmax,k,pk0,kn,sigma,W,skew,pktype)
  implicit none
  !I/O
  character(8), intent(in) :: cpmodel, model
  integer, intent(in) :: olmin, olmax, lmin, lmax, zn, kn
  double precision, intent(in) :: zs
  double precision, intent(in), dimension(1:2) :: sigma
  double precision, intent(in), dimension(1:zn) :: z, dz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(0:lmax) :: W
  double precision, intent(out), dimension(1:3,0:olmax) :: skew
  !optional
  character(4), intent(in), optional :: pktype
  !f2py character(4) :: pktype = 'T12'
  !internal
  character(4) :: pkfit = 'T12'
  integer :: oL(2), eL(2), tL(2)
  double precision, allocatable :: wp(:,:), ck(:,:)
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  ! cosmological parameters (to compute background quantities)
  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  tL(1) = 1
  tL(2) = max(lmax,olmax)
  oL = (/olmin,olmax/)
  eL = (/lmin,lmax/)

  ! skewspectrum
  allocate(wp(zn,tL(2)),ck(zn,tL(2)))
  if (present(pktype)) pkfit = pktype
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,tL,model,pkfit) 
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,tL,wp,ck)
  skew = 0d0
  write(*,*) 'calc skewspec'
  call skewspec_lens(cp,b,oL,eL,sigma,wp,ck,model,W(1:),skew(:,1:))
  deallocate(wp,ck)

end subroutine skewspeclens


end module bispec

