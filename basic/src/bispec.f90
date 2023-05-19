!///////////////////////////////////////////////////////////////////////!
! Lensing Bispectrum
!///////////////////////////////////////////////////////////////////////!

module bispec
  use general,   only: GLdxs, GLpoints, linspace, gauss_legendre_params, gl_initialize, gl_finalize
  use cosmofunc, only: cosmoparams, set_cosmoparams, omega_m
  use pstool,    only: binned_ells
  use bstool!, only: bispecfunc, bispec_lens_bin, bispec_lens_snr, bispec_lens_lss, bispec_lens_pb, bispec_lens_lss_init, bispec_lens_pb_init, zinterp, skewspec_lens, bispec_gauss_bin, cl_zweight, get_pl, cl_calc_flat, gal_zweight
  implicit none

contains


subroutine cl_flat(cpmodel,z,dz,zs,lmax,k,pk0,zn,kn,cl,pktype,cltype,dNdz,wdel)
!*  Compute flat-sky lensing powerspectrum analytically
!* 
!*  Args:
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs[2] (double) : two source redshifts
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :k[kn] (double) : k for the matter power spectrum in unit of [h/Mpc]
!*    :pk0[kn] (double) : the linear matter power spectrum at z=0 in unit of [Mpc^3/h^3]
!*
!*  Args(optional):
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*    :cltype (str) : kk, gk, or gg
!*    :dNdz[zn] (double) : redshift distribution of galaxy, only used when cltype includes g
!*    :wdel[zn,l] (double) : modified chi-kernel function for z-cleaning at l=0 to lmax
!*
!*  Returns:
!*    :cl[l] (double) : power spectrum from LSS contributions at [lmin,lmax]
!*
  implicit none
  !I/O
  character(8), intent(in) :: cpmodel, pktype, cltype
  integer, intent(in) :: lmax, zn, kn
  double precision, intent(in), dimension(2) :: zs
  double precision, intent(in), dimension(1:zn) :: z, dz, dNdz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(1:zn,0:lmax) :: wdel
  double precision, intent(out), dimension(0:lmax) :: cl
  !internal
  type(cosmoparams) :: cp
  double precision :: pl(zn,lmax), zker(zn), weight(2,zn,lmax)
  !opt4py :: zn = 0
  !opt4py :: kn = 0
  !add2py :: if zn == 0: zn=len(z)
  !add2py :: if kn == 0: kn=len(k)
  !opt4py :: pktype = 'T12'
  !opt4py :: cltype = 'kk'
  !opt4py :: dNdz = None
  !add2py :: if dNdz is None: dNdz = z*0.
  !opt4py :: wdel = None
  !add2py :: if wdel is None: wdel = numpy.zeros((zn,lmax+1))
  !add2py :: if len(dNdz) != zn: print('size of dNdz is strange')

  ! cosmological parameters (to compute background quantities)
  call set_cosmoparams(cp,cpmodel)

  ! precompute quantities for power spectrum
  call get_pl(cp,z,k*cp%h,pk0/cp%h**3,lmax,pktype,pl)
  call cl_zweight(cp,z,dz,zs,dNdz,wdel(:,1:),cltype,zker,weight)

  ! commpute power pectrum
  cl = 0d0
  call cl_calc_flat(pl,zker,weight,cl(1:))

end subroutine cl_flat


subroutine bispeclens(shap,cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,lan,kan,bl0,bl1,zn,kn,pktype,ltype,btype,dNdz,wdel)
!*  Compute lensing bispectrum analytically
!* 
!*  Args:
!*    :shap (str)  : shape of the bispectrum (equi, fold, sque, or isos)
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs[3] (double) : source redshifts
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :k[kn] (double) : k for the matter power spectrum in unit of [h/Mpc]
!*    :pk0[kn] (double) : the linear matter power spectrum at z=0 in unit of [Mpc^3/h^3]
!*
!*  Args(optional):
!*    :lan, kan (double) : parameters for the modified gravity extension, default to lan=kan=1 (GR)
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin=Linear, S02=Smith et al. 2002 or T12=Takahashi et al. 2012), default to T12
!*    :ltype (str) : fullsky correction (full) or not 
!*    :btype (str) : bispectrum type, i.e., kkk (lens-lens-lens), gkk (density-lens-lens), ggk (density-density-lens), or ggg (density-density-density)
!*    :dNdz[zn] (double) : redshift distribution of galaxy, only used when btype includes g
!*    :wdel[zn,l] (double) : modified chi-kernel function by z-cleaning at l=0 to lmax
!*
!*  Returns:
!*    :bl0[l] (double) : lensing bispectrum from LSS contributions at [lmin,lmax]
!*    :bl1[l] (double) : lensing bispectrum from post-Born contributions at [lmin,lmax]
!*
  implicit none
  !I/O
  character(8), intent(in) :: shap, cpmodel, model, pktype, ltype, btype
  integer, intent(in) :: lmin, lmax, zn, kn
  double precision, intent(in), dimension(3) :: zs
  double precision, intent(in) :: lan, kan
  double precision, intent(in), dimension(1:zn) :: z, dz, dNdz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(1:zn,0:lmax) :: wdel
  double precision, intent(out), dimension(lmin:lmax) :: bl0, bl1
  !internal
  integer :: l0, oL(2), i
  double precision, allocatable :: wp(:,:,:), wck(:,:,:,:), wdel_tmp(:,:)
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b
  !opt4py :: zn = 0
  !opt4py :: kn = 0
  !add2py :: if zn == 0: zn=len(z)
  !add2py :: if kn == 0: kn=len(k)
  !opt4py :: lan = 0.
  !opt4py :: kan = 0.
  !opt4py :: pktype = 'T12'
  !opt4py :: ltype = ''
  !opt4py :: btype = 'kkk'
  !opt4py :: dNdz = None
  !add2py :: if dNdz is None: dNdz = z*0.
  !opt4py :: wdel = None
  !add2py :: if wdel is None: wdel = numpy.zeros((zn,lmax+1))
  !add2py :: if len(dNdz) != zn: print('size of dNdz is strange')

  ! cosmological parameters (to compute background quantities)
  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  if (lmin<1) stop 'lmin should be >=1' 
  oL = (/lmin,lmax/)

  ! precompute quantities for bispectrum
  allocate(wp(3,zn,oL(2)),wck(3,3,zn,oL(2)),wdel_tmp(zn,1:lmax))
  wdel_tmp = wdel(:,1:)
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,oL,dNdz,wdel_tmp,model,pktype,btype) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,oL,b%weight,wp,wck) !only compute for btype=kkk
  deallocate(wdel_tmp)
  write(*,*) 'A'
  !if (btype/='kkk')  call gal_zweight(btype,b%zker,dNdz)

  ! MG parameter z-evolution
  allocate(b%mgp(2,zn));  b%mgp=1d0
  do i = 1, zn
    if (lan/=0d0) b%mgp(1,i) = (omega_m(1d0/(1d0+z(i)),cp))**lan
    if (kan/=0d0) b%mgp(2,i) = (omega_m(1d0/(1d0+z(i)),cp))**kan
  end do
  
  ! bispectrum (density and post-Bron bispectrum)
  write(*,*) 'compute bispectrum: shape = ', shap
  l0 = lmin+(lmax-lmin)*0.5d0/20d0
  if (mod(l0,2)/=0) l0 = l0+1
  call bispec_lens_lss(cp,b,shap,oL,model,ltype,l0,bl0)
  call bispec_lens_pb(shap,oL,wp,wck,bl1,l0,ltype)

  deallocate(b%mgp,wp,wck)

end subroutine bispeclens


subroutine bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan,bc,bl0,bl1,pktype,btype,dNdz,wdel)
!*  Compute binned lensing bispectrum analytically
!* 
!*  Args:
!*    :shap (str)  : shape of the bispectrum (equi, fold, sque, or isos)
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs[3] (double) : source redshifts
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :bn (int) : number of multipole bins
!*    :k[kn] (double) : k for the matter power spectrum in unit of [h/Mpc]
!*    :pk0[kn] (double) : the linear matter power spectrum at z=0 in unit of [Mpc^3/h^3]
!*
!*  Args(optional):
!*    :lan, kan (double) : parameters for the modified gravity extension, default to lan=kan=1 (GR)
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*    :btype (str) : bispectrum type, i.e., kkk (lens-lens-lens), gkk (density-lens-lens), ggk (density-density-lens), or ggg (density-density-density)
!*    :dNdz[zn] (double) : redshift distribution of galaxy, only used when btype includes g
!*    :wdel[zn,l] (double) : modified chi-kernel function by z-cleaning at l=0 to lmax
!*
!*  Returns:
!*    :bc[bn] (double)  : multipole bin centers
!*    :bl0[bn] (double) : binned lensing bispectrum from LSS contributions
!*    :bl1[bn] (double) : binned lensing bispectrum from post-Born contributions
!*
  implicit none
  !I/O
  character(8), intent(in) :: shap, cpmodel, model, pktype, btype
  integer, intent(in) :: lmin, lmax, zn, bn, kn
  double precision, intent(in), dimension(3) :: zs
  double precision, intent(in) :: lan, kan
  double precision, intent(in), dimension(1:zn) :: z, dz, dNdz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(1:zn,0:lmax) :: wdel
  double precision, intent(out), dimension(1:bn) :: bc, bl0, bl1
  !opt4py :: zn = 0
  !opt4py :: kn = 0
  !add2py :: if zn == 0: zn=len(z)
  !add2py :: if kn == 0: kn=len(k)
  !opt4py :: lan = 0.
  !opt4py :: kan = 0.
  !opt4py :: pktype = 'T12'
  !opt4py :: btype = 'kkk'
  !opt4py :: dNdz = None
  !add2py :: if dNdz is None: dNdz = z*0.
  !opt4py :: wdel = None
  !add2py :: if wdel is None: wdel = numpy.zeros((zn,lmax+1))
  !add2py :: if len(dNdz) != zn: print('size of dNdz is strange')

  !internal
  integer :: oL(2), l, i, eL1(2), eL2(2), eL3(2)
  double precision, allocatable :: bp(:), wp(:,:,:), wck(:,:,:,:)
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  if (lmin<1) stop 'lmin should be >=1' 
  oL = (/lmin,lmax/)

  ! precompute quantities for bispectrum
  allocate(wp(3,zn,oL(2)),wck(3,3,zn,oL(2)))
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,oL,dNdz,wdel(:,1:),model,pktype) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,oL,b%weight,wp,wck) !only compute for btype=kkk
  !if (btype/='kkk')  call gal_zweight(btype,b%zker,dNdz)

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
    case ('angl','isos')
      eL1 = int(bp(l:l+1))
      eL2 = int(bp(bn/2:bn/2+1))
      eL3 = int(bp(bn/2:bn/2+1))
    end select
    write(*,*) eL1, eL2, eL3
    call bispec_lens_bin(cp,b,eL1,eL2,eL3,wp,wck,model,bl0(l),bl1(l))
  end do

  deallocate(b%mgp,bp,wp,wck)

end subroutine bispeclens_bin


subroutine bispeclens_snr(cpmodel,model,z,dz,zn,zs,lmin,lmax,cl,k,pk0,kn,snr,pktype,btype,dNdz,cgg,ro,wdel)
!*  Compute SNR of lensing bispectrum analytically
!* 
!*  Args:
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, or modelp)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :z[zn] (double) : redshift points for the z-integral
!*    :dz[zn] (double) : interval of z
!*    :zs[3] (double) : source redshifts
!*    :lmin/lmax (int) : minimum/maximum multipoles of the bispectrum
!*    :cl[l] (int) : observed angular power spectrum at 0<=l<=lmax
!*    :k[kn] (double) : k for the matter power spectrum in unit of [h/Mpc]
!*    :pk0[kn] (double) : the linear matter power spectrum at z=0 in unit of [Mpc^3/h^3]
!*
!*  Args(optional):
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*    :btype (str) : bispectrum type, i.e., kkk (lens-lens-lens), gkk (density-lens-lens), ggk (density-density-lens), or ggg (density-density-density)
!*    :dNdz[zn] (double) : redshift distribution of galaxy, only used when btype includes g
!*    :wdel[zn,l] (double) : modified chi-kernel function by z-cleaning at l=0 to lmax
!*    :cgg[l] (double) : observed galaxy spectrum
!*    :ro (int) : output progress for every "ro" multipoles (ro=100, default)
!*
!*  Returns:
!*    :snr[2] (double) : total SNR amd LSS-only SNR
!*
  implicit none
  !I/O
  character(8), intent(in) :: cpmodel, model, pktype, btype
  integer, intent(in) :: lmin, lmax, zn, kn, ro
  double precision, intent(in), dimension(3) :: zs
  double precision, intent(in), dimension(1:zn) :: z, dz, dNdz
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(0:lmax) :: cl, cgg
  double precision, intent(in), dimension(1:zn,0:lmax) :: wdel
  double precision, intent(out), dimension(1:2) :: snr
  !opt4py :: zn = 0
  !opt4py :: kn = 0
  !add2py :: if zn == 0: zn=len(z)
  !add2py :: if kn == 0: kn=len(k)
  !opt4py :: pktype = 'T12'
  !opt4py :: btype = 'kkk'
  !opt4py :: dNdz = None
  !add2py :: if dNdz is None: dNdz = z*0.
  !add2py :: if len(dNdz) != zn: print('size of dNdz is strange')
  !opt4py :: wdel = None
  !add2py :: if wdel is None: wdel = numpy.zeros((zn,lmax+1))
  !opt4py :: cgg = None
  !opt4py :: ro = 100
  !add2py :: if cgg is None: cgg = cl*0.
  !add2py :: if len(cgg) != lmax+1: print('size of cgg is strange')

  !internal
  integer :: eL(2)
  double precision, allocatable :: wp(:,:,:), wck(:,:,:,:)
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  call set_cosmoparams(cp,cpmodel)

  ! other parameters
  if (lmin<1) stop 'lmin should be >=1' 
  eL = (/lmin,lmax/)

  ! precompute quantities for bispectrum
  allocate(wp(3,zn,eL(2)),wck(3,3,zn,eL(2)))
  call bispec_lens_lss_init(cp,b,z,dz,zs,k*cp%h,pk0/cp%h**3,eL,dNdz,wdel(:,1:),model,pktype,btype) !correction for h/Mpc to /Mpc
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zs,eL,b%weight,wp,wck) !only compute for btype=kkk
  !if (btype/='kkk')  call gal_zweight(btype,b%zker,dNdz)

  ! snr
  allocate(b%mgp(2,zn));  b%mgp=1d0
  write(*,*) 'compute bispectrum snr'
  select case(btype)
  case('kkk','ggg')
    call bispec_lens_snr(cp,b,eL,cl(1:lmax),wp,wck,model,ro,snr)
  case('gkk','ggk')
    call bispec_lens_snr_cross(cp,b,eL,cgg(1:lmax),cl(1:lmax),btype,model,ro,snr(1))
  end select
  deallocate(b%mgp,wp,wck)

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
    case ('angl','isos')
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


subroutine skewspeclens(cpmodel,model,zmin,zmax,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,skew,theta,pktype,pb,Om,H0,w0,wa,mnu,ns,verbose,wdel)
!* Compute skew spectrum using a matter bispectrum fitting formula
!*
!*  Args:
!*    :cpmodel (str) : cosmological parameter model (model0, modelw, modelp, or input)
!*    :model (str) : fitting formula of the matter bispectrum (LN=linear, SC=SC03, GM=Gil-Marin+12, 3B=3-shape-bispectrum, or RT=Takahashi+19)
!*    :zmin/zmax (double) : minimum/maximum z for z-integral
!*    :zn (int) : number of redshifts for the z-integral
!*    :zs[2] (double) : source redshifts where zs[2] is used for the squared map
!*    :lmin/lmax (int) : minimum/maximum multipoles of alms included in the skew spectrum
!*    :ols[bn] (int) : output multipoles to be computed for skew spectrum
!*    :k[kn] (double) : k for the matter power spectrum [h/Mpc]
!*    :pk0 (double) : the linear matter power spectrum at z=0 [Mpc^3/h^3]
!*
!*  Args(optional):
!*    :pktype (str) : fitting formula for the matter power spectrum (Lin, S02 or T12)
!*    :theta (double) : kappa map resolution in arcmin
!*    :pb (bool) : with post-Born correction or not (default=True)
!*    :verbose (bool) : output messages
!*    :wdel[zn,l] (double) : modified chi-kernel function by z-cleaning at l=0 to lmax
!*
!*  Returns:
!*    :skew (double)  : skew-spectrum
!*
  implicit none
  !I/O
  character(8), intent(in) :: cpmodel, model, pktype
  integer, intent(in) :: bn, lmin, lmax, zn, kn
  double precision, intent(in), dimension(2) :: zs
  double precision, intent(in) :: theta, zmin, zmax
  double precision, intent(in) :: Om, H0, w0, wa, mnu, ns
  integer, intent(in), dimension(1:bn) :: ols
  double precision, intent(in), dimension(1:kn) :: k, pk0
  double precision, intent(in), dimension(1:zn,0:lmax) :: wdel
  double precision, intent(out), dimension(1:3,1:bn) :: skew
  logical, intent(in) :: pb, verbose
  !opt4py :: pktype = 'T12'
  !opt4py :: pb = True
  !opt4py :: theta = 0.0
  !opt4py :: Om = 0.3
  !opt4py :: H0 = 70.
  !opt4py :: w0 = -1.
  !opt4py :: wa = 0.
  !opt4py :: mnu = 0.06
  !opt4py :: ns = 0.965
  !opt4py :: bn = 0
  !add2py :: if bn==0: bn=len(ols)
  !opt4py :: kn = 0
  !add2py :: if kn == 0: kn=len(k)
  !opt4py :: verbose = True
  !opt4py :: wdel = None
  !add2py :: if wdel is None: wdel = numpy.zeros((zn,lmax+1))
  !internal
  integer :: n, eL(2), tL(2)
  double precision :: zss(3)
  double precision, dimension(1:zn) :: z, dz, dndz
  double precision, allocatable :: wp(:,:,:), wck(:,:,:,:)
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp
  type(bispecfunc)  :: b

  ! cosmological parameters (to compute background quantities)
  if (cpmodel=='input') then
    cp%Om = Om
    cp%H0 = H0
    cp%w0 = w0
    cp%wa = wa
    cp%nu = mnu/(93.14d0*(cp%H0/100d0)**2*cp%Om)
    cp%h  = cp%H0/100d0
    cp%Ov = 1d0 - cp%Om
    cp%ns = ns
  else
    call set_cosmoparams(cp,cpmodel)
  end if

  ! z points for integral
  call zinterp(zmin,zmax,zn,1,z,dz)

  ! other parameters
  tL(1) = 1
  tL(2) = max(lmax,ols(bn))

  ! skewspectrum
  allocate(wp(3,zn,tL(2)),wck(3,3,zn,tL(2)))
  zss(1) = zs(1)
  zss(2) = zs(2)
  zss(3) = zs(2)
  dndz = 0d0 ! dummy
  call bispec_lens_lss_init(cp,b,z,dz,zss,k*cp%h,pk0/cp%h**3,tL,dndz,wdel(:,1:),model,pktype,verbose=verbose) 
  call bispec_lens_pb_init(cp,b%kl,b%pl,z,dz,zss,tL,b%weight,wp,wck)
  skew = 0d0
  if (verbose) write(*,*) 'calc skewspec'
  do n = 1, bn
    call skewspec_lens(cp,b,ols(n),(/lmin,lmax/),(/1d0,1d0/),wp,wck,model,theta,skew(:,n),pb)
    if (verbose) write(*,*) ols(n), skew(1,n)
  end do
  deallocate(wp,wck)

end subroutine skewspeclens


end module bispec

