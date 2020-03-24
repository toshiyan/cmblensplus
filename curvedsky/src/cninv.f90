!////////////////////////////////////////////////////!
! Map filtering
!////////////////////////////////////////////////////!

module cninv
  !from F90/src_utils
  use constants, only: pi
  use general, only: savetxt, check_error, str
  use hp_udgrade, only: udgrade_ring_1d_d
  use hp_cgd
  implicit none

  private pi

contains


subroutine cnfilter_freq(n,mn,npix,lmax,cl,bl,iNcov,maps,xlm,chn,lmaxs,nsides,itns,eps,fratio,filter,verbose,ro)
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
  integer, intent(in) :: n, mn, npix, lmax, chn, ro
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
  !opt4py :: ro = 50
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(mn,0:lmax) :: bl
  double precision, intent(in), dimension(n,mn,0:npix-1) :: iNcov, maps
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain)  :: mgc
  integer :: c, ni, mi, rn, nside, ilmaxs(chn), mnmaxs(chn), insides(chn,mn)
  integer(8) :: t1, t2, t_rate, t_max
  double precision :: clh(n,mn,0:lmax,0:lmax), ratio(itns(1))
  double precision, allocatable :: nij(:,:)
  double complex :: b(n,0:lmax,0:lmax)

  mnmaxs = mn

  !compute beam-convolved half signal spectrum
  call clhalf(n,mn,lmax,cl,bl,clh)

  !set multigrid parameters
  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides(1,:) = int(dsqrt(npix/12d0))
  else
    call check_error(lmax/=lmaxs(1),'input lmax is wrong',str(lmax)//','//str(lmaxs(1)))
    call check_error(npix/=12*nsides(1)**2,'input npix0 is wrong',str(npix)//','//str(12*nsides(1)**2))
    ilmaxs  = lmaxs
    do c = 1, chn
      insides(c,:) = nsides(c) 
    end do
    mnmaxs(2:) = 1
  end if
  call set_mgchain(mgc,chn,mn,mnmaxs,ilmaxs,insides,itns,eps,verbose,ro)

  !inverse noise covariance x map
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    do c = 1, mgc%n
      allocate(mgc%cv(c,mi)%nij(n,0:mgc%npix(c,mi)-1),mgc%cv(c,mi)%clh(n,0:mgc%lmax(c),0:mgc%lmax(c)))
      mgc%cv(c,mi)%nij = 0d0
      mgc%cv(c,mi)%clh = 0d0
    end do
    ! check Nij size for chain = 1
    call check_error(size(iNcov,dim=3)/=mgc%npix(1,mi),'iNcov size is wrong',str(size(iNcov,dim=3))//','//str(mgc%npix(1,mi)))
    mgc%cv(1,mi)%nij = iNcov(:,mi,:)*maps(:,mi,:)
  end do

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  !setup signal and noise covariance
  do mi = 1, mn
    mgc%cv(1,mi)%clh = clh(:,mi,:,:)
    mgc%cv(1,mi)%nij = iNcov(:,mi,:)
  end do

  !setup for multigrid preconditioner
  if (mgc%n>1) then
    write(*,*) 'degrade inv noise cov'
    do c = 2, mgc%n
      do mi = 1, mn
        do ni = 1, n
          call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,:),mgc%nside(1,mi),mgc%cv(c,mi)%nij(ni,:),mgc%nside(c,mi))
        end do
        if (mi==1) cycle
        call check_error(size(mgc%cv(c,mi)%nij,dim=2)/=mgc%npix(c,1),'iNcov size is inconsistent')
        mgc%cv(c,1)%nij = mgc%cv(c,1)%nij + mgc%cv(c,mi)%nij
      end do
      mgc%cv(c,1)%clh = sum(clh(:,:,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn)
    end do
    call coarse_invmatrix(n,mgc%lmax(mgc%n),mgc)
  end if

  !run kernel
  call system_clock(t1)
  call cg_algorithm(n,lmax,b,xlm,mgc,1,ratio)
  call system_clock(t2, t_rate, t_max) 
  write(*,*) "real time:", (t2-t1)/dble(t_rate)

  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,xlm)

end subroutine cnfilter_freq


subroutine cnfilter_freq_nside(n,mn0,mn1,npix0,npix1,lmax,cl,bl0,bl1,iNcov0,iNcov1,maps0,maps1,xlm,chn,lmaxs,nsides0,nsides1,itns,eps,fratio,filter,verbose,reducmn,ro)
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
  integer, intent(in) :: n, mn0, mn1, npix0, npix1, lmax, chn, ro, reducmn
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
  !opt4py :: reducmn = 0
  !opt4py :: ro = 50
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(mn0,0:lmax) :: bl0
  double precision, intent(in), dimension(mn1,0:lmax) :: bl1
  double precision, intent(in), dimension(n,mn0,0:npix0-1) :: iNcov0, maps0
  double precision, intent(in), dimension(n,mn1,0:npix1-1) :: iNcov1, maps1
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain)  :: mgc
  integer :: c, ni, mi, mn, rn, nside, ilmaxs(chn), mnmaxs(chn), insides(chn,mn0+mn1)
  integer(8) :: t1, t2, t_rate, t_max
  double precision :: clh(n,mn0+mn1,0:lmax,0:lmax), bl(mn0+mn1,0:lmax), ratio(itns(1))
  double precision, allocatable :: nij(:,:)
  double complex :: b(n,0:lmax,0:lmax)

  mn = mn0 + mn1
  mnmaxs = mn

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
    call check_error(npix1/=12*nsides1(1)**2,'input npix1 is wrong',str(npix1)//','//str(12*nsides1(1)**2))
    ilmaxs  = lmaxs
    do c = 1, chn
      insides(c,:mn0)   = nsides0(c) 
      insides(c,mn0+1:) = nsides1(c) 
    end do
    select case(reducmn)
    case (1)
      mnmaxs(2:) = 2
      insides(2:,2) = nsides1(2:)
    case (2)
      mnmaxs(2) = 2
      mnmaxs(3:) = 1
      insides(2,2) = nsides1(2)
    end select
  end if
  call set_mgchain(mgc,chn,mn,mnmaxs,ilmaxs,insides,itns,eps,verbose,ro)

  !inverse noise covariance x map
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    do c = 1, mgc%n
      allocate(mgc%cv(c,mi)%nij(n,0:mgc%npix(c,mi)-1),mgc%cv(c,mi)%clh(n,0:mgc%lmax(c),0:mgc%lmax(c)))
      mgc%cv(c,mi)%nij = 0d0
      mgc%cv(c,mi)%clh = 0d0
    end do
    ! check Nij size for chain = 1
    if (mi<=mn0) then
      call check_error(size(iNcov0,dim=3)/=mgc%npix(1,mi),'iNcov0 size is wrong',str(size(iNcov0,dim=3))//','//str(mgc%npix(1,mi)))
      mgc%cv(1,mi)%nij = iNcov0(:,mi,:)*maps0(:,mi,:)
    else
      call check_error(size(iNcov1,dim=3)/=mgc%npix(1,mi),'iNcov1 size is wrong',str(size(iNcov1,dim=3))//','//str(mgc%npix(1,mi)))
      mgc%cv(1,mi)%nij = iNcov1(:,mi-mn0,:)*maps1(:,mi-mn0,:)
    end if
  end do

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  !setup signal and noise covariance
  do mi = 1, mn
    mgc%cv(1,mi)%clh = clh(:,mi,:,:)
    if (mi<=mn0)  mgc%cv(1,mi)%nij = iNcov0(:,mi,:)
    if (mi>mn0)   mgc%cv(1,mi)%nij = iNcov1(:,mi-mn0,:)
  end do

  !setup for multigrid preconditioner
  if (mgc%n>1) then
    write(*,*) 'degrade inv noise cov'
    do c = 2, mgc%n

      select case(mgc%mnmax(c))
      case (4)
        do mi = 1, mn
          mgc%cv(c,mi)%clh = clh(:,mi,0:mgc%lmax(c),0:mgc%lmax(c))
          !store noise covariance for each sub-chain
          do ni = 1, n
            call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,:),mgc%nside(1,mi),mgc%cv(c,mi)%nij(ni,:),mgc%nside(c,mi))
          end do
        end do
      case (2)
        !for reduced number of mn
        do rn = 1, 2
            if (rn==1) nside = nsides0(c)
            if (rn==2) nside = nsides1(c)
            mgc%cv(c,rn)%nij = 0d0
            do mi = 1+(rn-1)*mn0, mn1 + (rn-1)*mn0
              allocate(nij(n,0:12*nside**2-1)); nij = 0d0
              do ni = 1, n
                call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,:),mgc%nside(1,mi),nij(ni,:),nside)
              end do
              mgc%cv(c,rn)%nij = mgc%cv(c,rn)%nij + nij
              deallocate(nij)
            end do
        end do
        mgc%cv(c,1)%clh = sum(clh(:,:mn0,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn0)
        mgc%cv(c,2)%clh = sum(clh(:,mn0+1:,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn1)
      case (1)
        do mi = 1, mn
          do ni = 1, n
            call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,:),mgc%nside(1,mi),mgc%cv(c,mi)%nij(ni,:),mgc%nside(c,mi))
          end do
        end do
        do mi = 2, mn
          if (size(mgc%cv(c,mi)%nij,dim=2)/=mgc%npix(c,1)) stop 'Nij size is inconsistent'
          mgc%cv(c,1)%nij = mgc%cv(c,1)%nij + mgc%cv(c,mi)%nij
        end do
        mgc%cv(c,1)%clh = sum(clh(:,:,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn)
      end select

    end do
    call coarse_invmatrix(n,mgc%lmax(mgc%n),mgc)
  end if

  !run kernel
  call system_clock(t1)
  call cg_algorithm(n,lmax,b,xlm,mgc,1,ratio)
  call system_clock(t2, t_rate, t_max) 
  write(*,*) "real time:", (t2-t1)/dble(t_rate)

  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,xlm)

end subroutine cnfilter_freq_nside


end module cninv

