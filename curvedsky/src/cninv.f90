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


subroutine cnfilter_freq(n,mn,npix,lmax,cl,bl,iNcov,maps,xlm,chn,lmaxs,nsides,itns,eps,filter,inl,verbose,ro,stat)
!* Combining multiple frequency maps optimally. 
!* The filtering would work if the noise variance is not significantly varied with scale (multipole). 
!*
!* Args:
!*    :n (int) : Number of maps, i.e., temperature only (n=1), polarization only (n=2) or both (n=3)
!*    :mn (int) : Number of frequencies
!*    :nside (int) : Nside of input map
!*    :lmax (int) : Maximum multipole of the input cl
!*    :cl[n,l] (double) : Theory signal power spectrum, with bounds (0:n-1,0:lmax)
!*    :bl[mn,l] (double) : Beam spectrum, with bounds (0:mn-1,0:lmax)
!*    :iNcov[n,mn,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:mn-1,0:npix-1)
!*    :maps[n,mn,pix] (double) : Input T, Q, U maps, with bouds (0:n-1,0:mn-1,0:npix-1)
!*
!* Args(optional):
!*    :chn (int) : number of grids for preconsitioner (chn=1 for diagonal preconditioner, default)
!*    :lmaxs[chain] (int) : Maximum multipole(s) at each preconditioning and lmaxs[0] is the input maximum multipole of cl
!*    :nsides[chain] (int) : Nside(s) of preconditoner and nsides[0] should be consistent with the input map's nside. 
!*    :eps[chain] (double): Numerical parameter to finish the iteration if ave(\|Ax-b\|)<eps, default to 1e-6
!*    :itns[chain] (int) : Number of interation(s)
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*    :inl[n,mn,l] (double) : inverse noise spectrum for non white noise case.
!*    :verbose (bool): Output messages, default to False
!*    :stat (str): Realtime status filename, default to no output file
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: stat
  integer, intent(in) :: n, mn, npix, lmax, chn, ro
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: verbose = False
  !opt4py :: ro = 50
  !opt4py :: stat = ''
  !opt4py :: inl = 0
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(mn,0:lmax) :: bl
  double precision, intent(in), dimension(n,mn,0:lmax) :: inl
  double precision, intent(in), dimension(n,mn,0:npix-1) :: iNcov, maps
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain)  :: mgc
  integer :: ou=6, c, ni, mi, rn, l, nside, ilmaxs(chn), mnmaxs(chn), insides(chn,mn)
  integer(8) :: t1, t2, t_rate, t_max
  double precision :: clh(n,mn,0:lmax,0:lmax)
  double precision, allocatable :: nij(:,:)
  double complex :: b(n,0:lmax,0:lmax)
  !replace
  !chargs :: npix -> nside
  !add2py :: npix = 12*nside**2

  mnmaxs = mn

  if (stat/='')  then 
    ou = 7
    open(unit=7,file=stat,status='replace')
  end if

  !compute beam-convolved half signal spectrum
  call clhalf(n,mn,lmax,cl,bl,clh)

  !set multigrid parameters
  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides(1,:) = int(dsqrt(npix/12d0))
  else
    call check_error(lmax/=lmaxs(1),'input lmax is wrong',str(lmax)//','//str(lmaxs(1)),ou=ou)
    call check_error(npix/=12*nsides(1)**2,'input npix0 is wrong',str(npix)//','//str(12*nsides(1)**2),ou=ou)
    ilmaxs  = lmaxs
    do c = 1, chn
      insides(c,:) = nsides(c) 
    end do
    mnmaxs(2:) = 1
  end if
  call set_mgchain(mgc,chn,mn,mnmaxs,ilmaxs,insides,itns,eps,verbose,ro)

  !initialize arrays for cinv
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    allocate( mgc%cv(1,mi)%imap(n,0:mgc%npix(1,mi)-1) )
    mgc%cv(1,mi)%imap = maps(:,mi,:)
    do c = 1, mgc%n
      allocate(mgc%cv(c,mi)%nij(n,0:mgc%npix(c,mi)-1),mgc%cv(c,mi)%clh(n,0:mgc%lmax(c),0:mgc%lmax(c)))
      mgc%cv(c,mi)%clh = 0d0
      mgc%cv(c,mi)%nij = 0d0
    end do
    ! check Nij size for chain = 1
    call check_error(size(iNcov,dim=3)/=mgc%npix(1,mi),'iNcov size is wrong',str(size(iNcov,dim=3))//','//str(mgc%npix(1,mi)),ou=ou)
    !setup signal and noise covariance
    mgc%cv(1,mi)%clh = clh(:,mi,:,:)
    mgc%cv(1,mi)%nij = iNcov(:,mi,:)
  end do

  !setup inverse noise spectrum if necessary
  if (sum(inl)/=0) then
    do mi = 1, mn
      do c = 1, mgc%n
        allocate(mgc%cv(c,mi)%nl(n,0:mgc%lmax(c),0:mgc%lmax(c)))
        mgc%cv(c,mi)%nl = 0d0
        do ni = 1, n
          do l = 2, mgc%lmax(c)
            mgc%cv(c,mi)%nl(ni,l,0:l)  = inl(ni,mi,l)
          end do
        end do
      end do
    end do
  end if

  !inverse noise covariance x map for computing b in rhs
  !do mi = 1, mn
    !mgc%cv(1,mi)%nij = iNcov(:,mi,:)*maps(:,mi,:)
  !end do

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)
  !call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  do mi = 1, mn
    deallocate(mgc%cv(1,mi)%imap)
  end do

  !setup for multigrid preconditioner
  if (mgc%n>1) then
    if (stat/=''.or.verbose)  write(ou,*) 'degrade inv noise cov'
    do c = 2, mgc%n
      do mi = 1, mn
        do ni = 1, n
          call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,:),mgc%nside(1,mi),mgc%cv(c,mi)%nij(ni,:),mgc%nside(c,mi))
        end do
        if (mi==1) cycle
        call check_error(size(mgc%cv(c,mi)%nij,dim=2)/=mgc%npix(c,1),'iNcov size is inconsistent',ou=ou)
        mgc%cv(c,1)%nij = mgc%cv(c,1)%nij + mgc%cv(c,mi)%nij
      end do
      mgc%cv(c,1)%clh = sum(clh(:,:,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn)
    end do
    call coarse_invmatrix(n,mgc%lmax(mgc%n),mgc)
  end if

  !run kernel
  call system_clock(t1)
  call cg_algorithm(n,lmax,b,xlm,mgc,1,ou)
  call system_clock(t2, t_rate, t_max) 

  if (stat/=''.or.verbose)  write(ou,*) "real time:", (t2-t1)/dble(t_rate)

  if (stat/='')  close(ou)

  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,xlm)

end subroutine cnfilter_freq


subroutine cnfilter_freq_nside(n,mn0,mn1,npix0,npix1,lmax,cl,bl0,bl1,iNcov0,iNcov1,maps0,maps1,xlm,chn,lmaxs,nsides0,nsides1,itns,eps,filter,inl,verbose,reducmn,ro,stat)
!* Same as cnfilter_freq but for the maps with two different Nsides. 
!*
!* Args:
!*    :n (int) : Number of maps, i.e., temperature only (n=1), polarization only (n=2) or both (n=3)
!*    :mn0/1 (int) : Number of frequencies
!*    :nside0/1 (int) : Nsides of input maps
!*    :lmax (int) : Maximum multipole of the input cl
!*    :cl[n,l] (double) : Theory signal power spectrum, with bounds (0:n-1,0:lmax)
!*    :bl0/1[mn,l] (double) : Beam function, with bounds (0:n-1,0:lmax)
!*    :iNcov0/1[n,mn,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :maps0/1[n,mn,pix] (double) : Input T, Q, U maps, with bouds (0:n-1,0:npix-1)
!*
!* Args(optional):
!*    :chn (int) : number of grids for preconsitioner (chn=1 for diagonal preconditioner, default)
!*    :lmaxs[chain] (int) : Maximum multipole(s) at each preconditioning and lmaxs[0] is the input maximum multipole of cl
!*    :nsides0/1[chain] (int) : Nside(s) of preconditoner and nsides[0] should be consistent with the input map's nside. 
!*    :eps[chain] (double): Numerical parameter to finish the iteration if ave(\|Ax-b\|)<eps, default to 1e-6
!*    :itns[chain] (int) : Number of interation(s)
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*    :verbose (bool): Output messages
!*    :reducmn (int): Reducing number of maps per chain (1,2) or not (0, default). If 1, the maps are combined for the same nside inside the multigrid chain. If 2, in addition to the procedure of 1, the each nside maprs are further combined into a single map inside the second chain (chain>=3).
!*    :stat (str): Realtime status filename
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: stat
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
  !opt4py :: verbose = False
  !opt4py :: reducmn = 0
  !opt4py :: ro = 50
  !opt4py :: stat = ''
  !opt4py :: inl = 0
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(mn0,0:lmax) :: bl0
  double precision, intent(in), dimension(mn1,0:lmax) :: bl1
  double precision, intent(in), dimension(n,mn0,0:npix0-1) :: iNcov0, maps0
  double precision, intent(in), dimension(n,mn1,0:npix1-1) :: iNcov1, maps1
  double precision, intent(in), dimension(n,mn0+mn1,0:lmax) :: inl
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mg_chain)  :: mgc
  integer :: ou=6, c, ni, mi, mn, rn, l, nside, ilmaxs(chn), mnmaxs(chn), insides(chn,mn0+mn1)
  integer(8) :: t1, t2, t_rate, t_max
  double precision :: clh(n,mn0+mn1,0:lmax,0:lmax), bl(mn0+mn1,0:lmax)
  double precision, allocatable :: nij(:,:)
  double complex :: b(n,0:lmax,0:lmax)

  integer :: nn

  !replace
  !chargs :: npix0 -> nside0
  !chargs :: npix1 -> nside1
  !add2py :: npix0 = 12*nside0**2
  !add2py :: npix1 = 12*nside1**2

  mn = mn0 + mn1
  mnmaxs = mn

  if (stat/='')  then 
    ou = 7
    open(unit=ou,file=stat,status='replace')
  end if

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
    call check_error(lmax/=lmaxs(1),'input lmax is wrong',str(lmax)//','//str(lmaxs(1)),ou)
    call check_error(npix0/=12*nsides0(1)**2,'input npix0 is wrong',str(npix0)//','//str(12*nsides0(1)**2),ou)
    call check_error(npix1/=12*nsides1(1)**2,'input npix1 is wrong',str(npix1)//','//str(12*nsides1(1)**2),ou)
    ilmaxs  = lmaxs
    ! set nsides for each freq map at each multigrid layer
    do c = 1, chn
      insides(c,:mn0)   = nsides0(c) 
      insides(c,mn0+1:) = nsides1(c) 
    end do
    ! change nsides for reduc freq and maps
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

  !setup signal and noise covariance
  allocate(mgc%cv(mgc%n,mn))
  do mi = 1, mn
    ! input map for rhs
    allocate( mgc%cv(1,mi)%imap(n,0:mgc%npix(1,mi)-1) )
    if (mi<=mn0)  mgc%cv(1,mi)%imap = maps0(:,mi,:)
    if (mi>mn0)   mgc%cv(1,mi)%imap = maps1(:,mi,:)
    ! signal and noise covariance
    do c = 1, mgc%n
      allocate( mgc%cv(c,mi)%nij(n,0:mgc%npix(c,mi)-1), mgc%cv(c,mi)%clh(n,0:mgc%lmax(c),0:mgc%lmax(c)) )
      mgc%cv(c,mi)%nij = 0d0
      mgc%cv(c,mi)%clh = 0d0
    end do
    ! for chain=1
    mgc%cv(1,mi)%clh = clh(:,mi,:,:)
    ! check Nij size for chain = 1
    if (mi<=mn0) then
      call check_error(size(iNcov0,dim=3)/=mgc%npix(1,mi),'iNcov0 size is wrong',str(size(iNcov0,dim=3))//','//str(mgc%npix(1,mi)),ou)
      mgc%cv(1,mi)%nij = iNcov0(:,mi,:)!*maps0(:,mi,:)
    else
      call check_error(size(iNcov1,dim=3)/=mgc%npix(1,mi),'iNcov1 size is wrong',str(size(iNcov1,dim=3))//','//str(mgc%npix(1,mi)),ou)
      mgc%cv(1,mi)%nij = iNcov1(:,mi-mn0,:)!*maps1(:,mi-mn0,:)
    end if
  end do

  !setup inverse noise spectrum if necessary
  if (sum(inl)/=0) then
    do mi = 1, mn
      do c = 1, mgc%n
        allocate(mgc%cv(c,mi)%nl(n,0:mgc%lmax(c),0:mgc%lmax(c)))
        mgc%cv(c,mi)%nl = 0d0
        do ni = 1, n
          do l = 2, mgc%lmax(c)
            mgc%cv(c,mi)%nl(ni,l,0:l)  = inl(ni,mi,l)
          end do
        end do
      end do
    end do
  end if

  !first compute b = C^1/2 N^-1 X
  call matmul_rhs(n,mn,mgc%npix(1,:),lmax,clh,mgc%cv(1,:),b)

  do mi = 1, mn
    deallocate(mgc%cv(1,mi)%imap)
  end do

  !do mi = 1, mn
  !  if (mi<=mn0)  mgc%cv(1,mi)%nij = iNcov0(:,mi,:)
  !  if (mi>mn0)   mgc%cv(1,mi)%nij = iNcov1(:,mi-mn0,:)
  !end do

  !setup for multigrid preconditioner
  ! for chain>=2
  if (mgc%n>1) then

    if (stat/='' .or. verbose) write(ou,*) 'degrade inv noise cov'
    
    do c = 2, mgc%n

      select case(mgc%mnmax(c))
      
      case (4:)

        if (mgc%mnmax(c)/=mn)  stop 'Inconsistency in the number of maps in the multigrid chain.'
      
        do mi = 1, mn
          mgc%cv(c,mi)%clh = clh(:,mi,0:mgc%lmax(c),0:mgc%lmax(c))
          !store noise covariance for each sub-chain
          do ni = 1, n
            call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,:),mgc%nside(1,mi),mgc%cv(c,mi)%nij(ni,:),mgc%nside(c,mi))
          end do
        end do
      
      case (2)  !for the case combining frequency maps

        do rn = 1, 2
            if (rn==1) nside = nsides0(c)
            if (rn==2) nside = nsides1(c)
            mgc%cv(c,rn)%nij = 0d0
            do mi = 1+(rn-1)*mn0, mn0+(rn-1)*mn1 !lopp for freqs
              allocate(nij(n,0:12*nside**2-1)); nij = 0d0
              do ni = 1, n !loop for tqu
                call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,0:),mgc%nside(1,mi),nij(ni,0:),nside)
              end do
              mgc%cv(c,rn)%nij = mgc%cv(c,rn)%nij + nij
              deallocate(nij)
            end do
        end do
        mgc%cv(c,1)%clh = sum(clh(:,:mn0,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn0)
        mgc%cv(c,2)%clh = sum(clh(:,mn0+1:,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn1)

      case (1)  ! for the case combining frequency and Nside maps to a single map
        
        ! downgrade resolution from the original to inside chain
        nside = mgc%nside(c,1)

        ! sum inverse noise covariance for freqs and Nsides
        mgc%cv(c,1)%nij = 0d0
        do mi = 1, mn !lopp for freqs
           allocate(nij(n,0:12*nside**2-1)); nij = 0d0
           do ni = 1, n !loop for tqu
              call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,0:),mgc%nside(1,mi),nij(ni,0:),nside)
           end do
           mgc%cv(c,1)%nij = mgc%cv(c,1)%nij + nij
           deallocate(nij)
        end do
        !do mi = 1, mn
        !  do ni = 1, n
        !    call udgrade_ring_1d_d(mgc%cv(1,mi)%nij(ni,0:),mgc%nside(1,mi),mgc%cv(c,mi)%nij(ni,0:),mgc%nside(c,mi))
        !  end do
        !end do
        !do mi = 2, mn
        !  call check_error(size(mgc%cv(c,mi)%nij,dim=2)/=mgc%npix(c,mi),'Nij is inconsistent',str(size(mgc%cv(c,mi)%nij,dim=2))//','//str(mgc%npix(c,mi)),ou)
        !  mgc%cv(c,1)%nij = mgc%cv(c,1)%nij + mgc%cv(c,mi)%nij
        !end do

        mgc%cv(c,1)%clh = sum(clh(:,:,0:mgc%lmax(c),0:mgc%lmax(c)),dim=2)/dble(mn)
      
      case default

        stop 'Specified multigrid chain is not supported. Please change the reducmn.'

      end select

    end do
    
    call coarse_invmatrix(n,mgc%lmax(mgc%n),mgc)
  
  end if

  !run kernel
  call system_clock(t1)
  call cg_algorithm(n,lmax,b,xlm,mgc,1,ou)
  call system_clock(t2, t_rate, t_max) 
  
  if (stat/='' .or. verbose)  write(ou,*) "real time:", (t2-t1)/dble(t_rate)

  if (stat/='')  close(ou)

  call free_mgchain(mgc)

  call correct_filtering(n,lmax,cl,filter,xlm)

end subroutine cnfilter_freq_nside


end module cninv

