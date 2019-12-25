!///////////////////////////////////////////////////////////////////////!
! Lensing Bispectrum
!///////////////////////////////////////////////////////////////////////!

module bispec
  use general,   only: GLdxs, GLpoints, linspace, gauss_legendre_params, gl_initialize, gl_finalize
  use cosmofunc, only: cosmoparams, set_cosmoparams, omega_m
  use pstool,    only: binned_ells
  use bstool,    only: bispecfunc, bispec_lens_bin, bispec_lens_snr, bispec_lens_lss, bispec_lens_pb, bispec_lens_lss_init, bispec_lens_pb_init, zinterp, skewspec_lens, bispec_gauss_bin
  implicit none

contains


subroutine bispeclens(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan,bl0,bl1,pktype,ltype)
!*  Compute lensing bispectrum analytically
!* 
!*  Args:
!*    :shap (str)  : shape of the bispectrum (equi, fold, sque, or isos)
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :zn (int) : number of redshifts for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs (double) : source redshift
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :k[kn] (double) : k for the matter power spectrum
!*    :pk0 (double) : the linear matter power spectrum at z=0
!*    :kn (int) : size of k
!*
!*  Args(optional):
!*    :lan, kan (double) : parameters for the modified gravity extension, default to lan=kan=1 (GR)
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*    :ltype (str) : fullsky correction (full) or not 
!*
!*  Returns:
!*    :bl0[l] (double) : lensing bispectrum from LSS contributions at [lmin,lmax]
!*    :bl1[l] (double) : lensing bispectrum from post-Born contributions at [lmin,lmax]
!*
  implicit none
  !I/O
  character(8), intent(in) :: shap, cpmodel, model, pktype, ltype
  integer, intent(in) :: lmin, lmax, zn, kn
  double precision, intent(in) :: lan, kan, zs
  double precision, intent(in), dimension(1:zn) :: z, dz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(out), dimension(lmin:lmax) :: bl0, bl1
  !opt4py :: lan = 0.
  !opt4py :: kan = 0.
  !opt4py :: pktype = 'T12'
  !opt4py :: ltype = ''
  !internal
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
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,oL,model,pktype) !correction for h/Mpc to /Mpc
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
  call bispec_lens_lss(cp,b,shap,oL,model,ltype,l0,bl0)
  call bispec_lens_pb(shap,oL,wp,ck,bl1,l0,ltype)

end subroutine bispeclens


subroutine bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,bc,bl0,bl1,pktype)
!*  Compute binned lensing bispectrum analytically
!* 
!*  Args:
!*    :shap (str)  : shape of the bispectrum (equi, fold, sque, or isos)
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :zn (int) : number of redshifts for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs (double) : source redshift
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :bn (int) : number of multipole bins
!*    :k[kn] (double) : k for the matter power spectrum
!*    :pk0 (double) : the linear matter power spectrum at z=0
!*    :kn (int) : size of k
!*
!*  Args(optional):
!*    :lan, kan (double) : parameters for the modified gravity extension, default to lan=kan=1 (GR)
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*
!*  Returns:
!*    :bc[bn] (double)  : multipole bin centers
!*    :bl0[bn] (double) : binned lensing bispectrum from LSS contributions
!*    :bl1[bn] (double) : binned lensing bispectrum from post-Born contributions
!*
  implicit none
  !I/O
  character(8), intent(in) :: shap, cpmodel, model, pktype
  integer, intent(in) :: lmin, lmax, zn, bn, kn
  double precision, intent(in) :: lan, kan, zs
  double precision, intent(in),  dimension(1:zn) :: z, dz
  double precision, intent(in),  dimension(1:kn) :: k, pk0
  double precision, intent(out), dimension(1:bn) :: bc, bl0, bl1
  !opt4py :: lan = 0.
  !opt4py :: kan = 0.
  !opt4py :: pktype = 'T12'
  !internal
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
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,oL,model,pktype) !correction for h/Mpc to /Mpc
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


subroutine bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,snr,pktype)
!*  Compute SNR of lensing bispectrum analytically
!* 
!*  Args:
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :zn (int) : number of redshifts for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs (double) : source redshift
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :cl[l] (int) : observed lensing spectrum at 0<=l<=lmax
!*    :k[kn] (double) : k for the matter power spectrum
!*    :pk0 (double) : the linear matter power spectrum at z=0
!*    :kn (int) : size of k
!*
!*  Args(optional):
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*
!*  Returns:
!*    :snr (double)  : total SNR
!*
  implicit none
  !I/O
  character(8), intent(in) :: cpmodel, model, pktype
  integer, intent(in) :: lmin, lmax, zn, kn
  double precision, intent(in) :: zs
  double precision, intent(in), dimension(1:zn) :: z, dz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(0:lmax) :: cl
  double precision, intent(out) :: snr
  !opt4py :: pktype = 'T12'
  !internal
  integer :: eL(2)
  double precision, allocatable :: wp(:,:), ck(:,:)
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  if (lmin<1) stop 'lmin should be >=1' 
  eL = (/lmin,lmax/)

  ! precompute quantities for bispectrum
  allocate(wp(zn,eL(2)),ck(zn,eL(2)))
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,eL,model,pktype) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,eL,wp,ck)

  ! snr
  allocate(b%mgp(2,zn));  b%mgp=1d0
  write(*,*) 'compute bispectrum snr'
  call bispec_lens_snr(cp,b,eL,cl(1:lmax),wp,ck,model,snr)
  deallocate(b%mgp,wp,ck)

end subroutine bispeclens_snr


subroutine bispeclens_gauss_bin(shap,bn,lmin,lmax,cl,bc,bl)
!*  Compute binned bispectrum analytically for the quadratic gaussian model
!*
!*  Args:
!*    :shap (str)      : shape of the bispectrum (equi, fold, sque, or isos)
!*    :bn (int)        : number of multipole bins
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :cl[l] (double)  : the power spectrum at [0:lmax+1]
!*
!*  Returns:
!*    :bc[bn] (double) : multipole bin centers
!*    :bl[bn] (double) : binned bispectrum
!*
  implicit none
  !I/O
  character(4), intent(in) :: shap
  integer, intent(in) :: lmin, lmax, bn
  double precision, intent(in), dimension(lmax+1) :: cl
  double precision, intent(out), dimension(bn) :: bc, bl
  !internal
  integer :: oL(2), l, i, eL1(2), eL2(2), eL3(2)
  double precision :: bp(bn+1)

  if (lmin<1) stop 'lmin should be >=1' 
  oL = (/lmin,lmax/)

  ! binning
  call binned_ells(oL,bp,bc)

  ! bispectrum
  write(*,*) 'compute bispectrum', shap
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

    call bispec_gauss_bin(eL1,eL2,eL3,cl,bl(bn))

  end do


end subroutine bispeclens_gauss_bin


subroutine zpoints(zmin,zmax,zn,z,dz,zspace)
!*  Precomputing interpolation points for z
!*
!*  Args: 
!*    :zmin/zmax (double) : minimum/maximum redshifts
!*    :zn (int)           : number of redshifts
!*
!*  Args(optional):
!*    :zspace (int) : type of spacing
!*
!*  Returns:
!*    :z[zn] (double)  : redshifts
!*    :dz[zn] (double) : redshift intervals
!*
  implicit none
  !I/O
  integer, intent(in) :: zn, zspace
  double precision, intent(in) :: zmin, zmax
  double precision, intent(out), dimension(1:zn) :: z, dz
  !opt4py :: zspace = 1

  call zinterp(zmin,zmax,zn,zspace,z,dz)

end subroutine zpoints


subroutine skewspeclens(cpmodel,model,z,dz,zn,zs,olmin,olmax,lmin,lmax,k,pk0,kn,sigma,W,skew,pktype)
!*  Compute skewspectrum analytically
  implicit none
  !I/O
  character(8), intent(in) :: cpmodel, model, pktype
  integer, intent(in) :: olmin, olmax, lmin, lmax, zn, kn
  double precision, intent(in) :: zs
  double precision, intent(in), dimension(1:2) :: sigma
  double precision, intent(in), dimension(1:zn) :: z, dz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(0:lmax) :: W
  double precision, intent(out), dimension(1:3,0:olmax) :: skew
  !opt4py :: pktype = 'T12'
  !internal
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
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,tL,model,pktype) 
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,tL,wp,ck)
  skew = 0d0
  write(*,*) 'calc skewspec'
  call skewspec_lens(cp,b,oL,eL,sigma,wp,ck,model,W(1:),skew(:,1:))
  deallocate(wp,ck)

end subroutine skewspeclens


end module bispec

