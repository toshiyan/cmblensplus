!////////////////////////////////////////////////////!
! Map filtering
!////////////////////////////////////////////////////!

module cninv
  !from F90/src_utils
  use constants, only: pi
  use general, only: savetxt, check_error, str
  use hp_cgd
  implicit none

  private pi

contains


subroutine cnfilter(n,npix,lmax,cl,bl,iNcov,maps,xlm,chn,lmaxs,nsides,itns,eps,fratio,filter,verbose)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is S+N
!*
!* Args:
!*    :n (int) : Number of maps, i.e., temperature only (n=1), polarization only (n=2) or both (n=3)
!*    :npix (int) : Number of pixels of input map(s)
!*    :lmax (int) : Maximum multipole of the input cl
!*    :cl[n,l] (double) : Theory signal power spectrum, with bounds (0:n-1,0:lmax)
!*    :bl[l] (double) : Beam spectrum with bounds (0:lmax)
!*    :iNcov[n,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :maps[n,pix] (double) : Input T, Q, U maps, with bouds (0:n-1,0:npix-1)
!*
!* Args(optional):
!*    :chn (int) : number of grids for preconsitioner (chn=1 for diagonal preconditioner, default)
!*    :lmaxs[chain] (int) : Maximum multipole(s) at each preconditioning and lmaxs[0] is the input maximum multipole of cl
!*    :nsides[chain] (int) : Nside(s) of preconditoner and nsides[0] should be consistent with the input map's nside. 
!*    :eps[chain] (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :itns[chain] (int) : Number of interation(s)
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*    :fratio (str): Output filename of |r|^2/|b^2|
!*    :verbose (bool): Check the matrix at the coarsest grid. 
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: fratio
  integer, intent(in) :: n, npix, lmax, chn
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: fratio = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(0:lmax) :: bl
  double precision, intent(in), dimension(n,0:npix-1) :: iNcov, maps
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain) :: mgc
  integer :: c, mi, mn, ilmaxs(chn), insides(chn,1)
  integer :: t1, t2, t_rate, t_max
  double precision, dimension(n,1,0:lmax,0:lmax) :: clh
  double precision, dimension(1,0:lmax) :: ibl
  double complex :: b(n,0:lmax,0:lmax), alm(n,0:lmax,0:lmax)
  double precision :: ratio(itns(1))

  mn = 1

  !compute beam-convolved half signal spectrum
  ibl(1,:) = bl
  call clhalf(n,mn,lmax,cl,ibl,clh)

  !set multigrid parameters
  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides = reshape((/int(dsqrt(npix/12d0))/),(/chn,1/))
  else
    call check_error(lmax/=lmaxs(1),'input lmax is wrong',str(lmax)//','//str(lmaxs(1)))
    call check_error(npix/=12*nsides(1)**2,'input npix is wrong',str(npix)//','//str(12*nsides(1)**2))
    ilmaxs  = lmaxs
    insides = reshape(nsides,(/chn,1/))
  end if
  call set_mgchain(mgc,chn,mn,ilmaxs,insides,itns,eps,verbose)

  !inverse noise covariance
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    do c = 1, mgc%n
      allocate(mgc%cv(c,mi)%nij(n,mgc%npix(c,mi)))
    end do
    if (size(iNcov,dim=2)/=mgc%npix(1,mi)) stop 'incov size is strange'
    mgc%cv(1,mi)%nij = iNcov*maps
  end do

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  do mi = 1, mn
    mgc%cv(1,mi)%nij = iNcov
  end do

  call system_clock(t1)
  call cg_algorithm(n,mn,lmax,clh,b,alm,mgc,1,ratio)
  call system_clock(t2, t_rate, t_max) 
  write(*,*) "real time:", (t2-t1)/dble(t_rate)
  
  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,alm)
  xlm = alm

end subroutine cnfilter


subroutine cnfilter_freq(n,mn,npix,lmax,cl,bl,iNcov,maps,xlm,chn,lmaxs,nsides,itns,eps,fratio,filter,verbose)
!* Same as cnfilter but combining multiple frequency maps which have a same nside. 
!*
!* Args:
!*    :n (int) : Number of maps, i.e., temperature only (n=1), polarization only (n=2) or both (n=3)
!*    :mn (int) : Number of frequencies
!*    :npix (int) : Number of pixels of input map(s)
!*    :lmax (int) : Maximum multipole of the input cl
!*    :cl[n,l] (double) : Theory signal power spectrum, with bounds (0:n-1,0:lmax)
!*    :bl[mn,l] (double) : Beam spectrum, with bounds (0:n-1,0:lmax)
!*    :iNcov[n,mn,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :maps[n,mn,pix] (double) : Input T, Q, U maps, with bouds (0:n-1,0:npix-1)
!*
!* Args(optional):
!*    :chn (int) : number of grids for preconsitioner (chn=1 for diagonal preconditioner, default)
!*    :lmaxs[chain] (int) : Maximum multipole(s) at each preconditioning and lmaxs[0] is the input maximum multipole of cl
!*    :nsides[chain] (int) : Nside(s) of preconditoner and nsides[0] should be consistent with the input map's nside. 
!*    :eps[chain] (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :itns[chain] (int) : Number of interation(s)
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*    :fratio (str): Output filename of |r|^2/|b^2|
!*    :verbose (bool): Check the matrix at the coarsest grid. 
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: fratio
  integer, intent(in) :: n, mn, npix, lmax, chn
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: fratio = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(mn,0:lmax) :: bl
  double precision, intent(in), dimension(n,mn,0:npix-1) :: iNcov, maps
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain) :: mgc
  integer :: c, mi, ilmaxs(chn), insides(chn,mn)
  integer(8) :: t1, t2, t_rate, t_max
  double precision, dimension(n,mn,0:lmax,0:lmax) :: clh
  double complex :: b(n,0:lmax,0:lmax)
  double precision :: ratio(itns(1))

  !compute beam-convolved half signal spectrum
  call clhalf(n,mn,lmax,cl,bl,clh)

  !set multigrid parameters
  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides(1,:) = int(dsqrt(npix/12d0))
  else
    call check_error(lmax/=lmaxs(1),'input lmax is wrong',str(lmax)//','//str(lmaxs(1)))
    call check_error(npix/=12*nsides(1)**2,'input npix is wrong',str(npix)//','//str(12*nsides(1)**2))
    ilmaxs  = lmaxs
    do c = 1, chn
      insides(c,:) = nsides(c) !use same nside for all frequencies
    end do
  end if
  call set_mgchain(mgc,chn,mn,ilmaxs,insides,itns,eps,verbose)

  !inverse noise covariance
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    do c = 1, mgc%n
      allocate(mgc%cv(c,mi)%nij(n,mgc%npix(c,mi)))
    end do
    call check_error(size(iNcov,dim=3)/=mgc%npix(1,mi),'iNcov size is strange',str(size(iNcov,dim=3))//','//str(mgc%npix(1,mi)))
    mgc%cv(1,mi)%nij = iNcov(:,mi,:)*maps(:,mi,:)
  end do

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  do mi = 1, mn
    mgc%cv(1,mi)%nij = iNcov(:,mi,:)
  end do

  call system_clock(t1)
  call cg_algorithm(n,mn,lmax,clh,b,xlm,mgc,1,ratio)
  call system_clock(t2, t_rate, t_max) 
  write(*,*) "real time:", (t2-t1)/dble(t_rate)

  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,xlm)

end subroutine cnfilter_freq


subroutine cnfilter_freq_nside(n,mn0,mn1,npix0,npix1,lmax,cl,bl0,bl1,iNcov0,iNcov1,maps0,maps1,xlm,chn,lmaxs,nsides0,nsides1,itns,eps,fratio,filter,verbose)
!* Same as cnfilter but combining multiple frequency maps and these maps are divided into two different nside groups. 
!*
!* Args:
!*    :n (int) : Number of maps, i.e., temperature only (n=1), polarization only (n=2) or both (n=3)
!*    :mn0/1 (int) : Number of frequencies
!*    :npix0/1 (int) : Number of pixels of input map(s)
!*    :lmax (int) : Maximum multipole of the input cl
!*    :cl[n,l] (double) : Theory signal power spectrum, with bounds (0:n-1,0:lmax)
!*    :bl0/1[mn,l] (double) : Beam spectrum, with bounds (0:n-1,0:lmax)
!*    :iNcov0/1[n,mn,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :maps0/1[n,mn,pix] (double) : Input T, Q, U maps, with bouds (0:n-1,0:npix-1)
!*
!* Args(optional):
!*    :chn (int) : number of grids for preconsitioner (chn=1 for diagonal preconditioner, default)
!*    :lmaxs[chain] (int) : Maximum multipole(s) at each preconditioning and lmaxs[0] is the input maximum multipole of cl
!*    :nsides0/1[chain] (int) : Nside(s) of preconditoner and nsides[0] should be consistent with the input map's nside. 
!*    :eps[chain] (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :itns[chain] (int) : Number of interation(s)
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*    :fratio (str): Output filename of |r|^2/|b^2|
!*    :verbose (bool): Check the matrix at the coarsest grid. 
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: fratio
  integer, intent(in) :: n, mn0, mn1, npix0, npix1, lmax, chn
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides0, nsides1, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides0 = [0]
  !opt4py :: nsides1 = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: fratio = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(mn0,0:lmax) :: bl0
  double precision, intent(in), dimension(mn1,0:lmax) :: bl1
  double precision, intent(in), dimension(n,mn0,0:npix0-1) :: iNcov0, maps0
  double precision, intent(in), dimension(n,mn1,0:npix1-1) :: iNcov1, maps1
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain) :: mgc
  integer :: c, mi, mn, ilmaxs(chn), insides(chn,mn0+mn1)
  integer :: t1, t2, t_rate, t_max
  double precision :: clh(n,mn0+mn1,0:lmax,0:lmax), bl(mn0+mn1,0:lmax)
  double precision :: ratio(itns(1))
  double complex :: b(n,0:lmax,0:lmax)

  mn = mn0 + mn1

  !compute beam-convolved half signal spectrum
  bl(:mn0,:)   = bl0
  bl(mn0+1:,:) = bl1
  call clhalf(n,mn,lmax,cl,bl,clh)

  !set multigrid parameters
  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides(1,:mn0)   = int(dsqrt(npix0/12d0))
    insides(1,mn0+1:) = int(dsqrt(npix1/12d0))
  else
    call check_error(lmax/=lmaxs(1),'input lmax is wrong',str(lmax)//','//str(lmaxs(1)))
    call check_error(npix0/=12*nsides0(1)**2,'input npix0 is wrong',str(npix0)//','//str(12*nsides0(1)**2))
    call check_error(npix1/=12*nsides1(1)**2,'input npix0 is wrong',str(npix1)//','//str(12*nsides1(1)**2))
    ilmaxs  = lmaxs
    do c = 1, chn
      insides(c,:mn0)   = nsides0(c) 
      insides(c,mn0+1:) = nsides1(c) 
    end do
  end if
  call set_mgchain(mgc,chn,mn,ilmaxs,insides,itns,eps,verbose)

  !inverse noise covariance
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    do c = 1, mgc%n
      allocate(mgc%cv(c,mi)%nij(n,mgc%npix(c,mi)))
    end do
    if (mi<=mn0) then
      if (size(iNcov0,dim=3)/=mgc%npix(1,mi)) stop 'iNcov size is strange'
      mgc%cv(1,mi)%nij = iNcov0(:,mi,:)*maps0(:,mi,:)
    else
      if (size(iNcov1,dim=3)/=mgc%npix(1,mi)) stop 'iNcov size is strange'
      mgc%cv(1,mi)%nij = iNcov1(:,mi-mn0,:)*maps1(:,mi-mn0,:)
    end if
  end do

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  do mi = 1, mn
    if (mi<=mn0) then
      mgc%cv(1,mi)%nij = iNcov0(:,mi,:)
    else
      mgc%cv(1,mi)%nij = iNcov1(:,mi-mn0,:)
    end if
  end do

  call system_clock(t1)
  call cg_algorithm(n,mn,lmax,clh,b,xlm,mgc,1,ratio)
  call system_clock(t2, t_rate, t_max) 
  write(*,*) "real time:", (t2-t1)/dble(t_rate)

  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,xlm)

end subroutine cnfilter_freq_nside


end module cninv

