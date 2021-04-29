! Simple program demonstrating how to generate a simulated lensed map
! AL, Feb 2004; Updated Oct 2007

program SimLensCMB

  use HealpixObj
  use HealpixVis
  use Random
  use spinalm_tools
  use IniFile
  use AMLUtils

  !namikawa
  use general,   only: str, savetxt, linspace
  use pstool,     only: calccl
  use healpix_types
  use alm_tools
  use pix_tools
  use fitstools
  use head_fits, only: add_card, write_minimal_header
  !

  implicit none
  Type(HealpixInfo)  :: H
  Type(HealpixMap)   :: M, GradPhi
  Type(HealpixPower) :: P
  Type(HealpixAlm)   :: A
  integer            :: nside, lmax, npix
  character(LEN=256)  :: w8name = '../Healpix_2.00/data/'
  character(LEN=256)  :: file_stem, f1, f2, f3, out_file_root, cls_lensed_file
  integer, parameter :: lens_interp =1, lens_exact = 2
  integer :: lens_method = lens_interp
  integer :: mpi_division_method = division_equalrows
  integer ::  interp_method,  rand_seed
  logical :: err, want_pol
  real :: interp_factor

  !namikawa 
  character(LEN=64) InLine
  logical :: WantGMAP, WantLCl, lensfix, TBEB
  integer :: n, i, ii, sn(2), simn
  double precision, parameter :: Tcmb = 2.72d6, pi0 = 3.14d0
  double precision, allocatable :: mCl(:,:), cl(:,:,:)
  complex(SPC), allocatable :: Grad(:,:)
  !

#ifdef MPIPIX
  call mpi_init(i)
#endif

  Ini_Fail_On_Not_Found = .true.
  call Ini_Open(GetParam(1), 3,err)
  if (err) then
#ifdef MPIPIX
    call mpi_finalize(i)
#endif
    stop 'No ini'
  end if
  nside  = Ini_Read_Int('nside')
  npix = 12*nside**2
  lmax   = Ini_Read_Int('lmax')  
  WantLCl = Ini_Read_Logical('WantLCl')
  f1 = Ini_Read_String('cls_file')
  TBEB = Ini_Read_Logical('TBEB') 
  out_file_root = Ini_Read_String('out_file_root')
  lens_method = Ini_Read_Int('lens_method')
  want_pol = Ini_Read_Logical('want_pol')
  rand_seed = Ini_Read_Int('rand_seed')
  interp_method = Ini_read_int('interp_method')
  Ini_Fail_On_Not_Found = .false.
  w8name = Ini_Read_String('w8dir')
  interp_factor=0
  if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
  mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif 

  !//// namikawa ////!
  WantGMAP = Ini_Read_Logical('WantGMAP')
  lensfix = Ini_Read_Logical('lensfix')
  sn(1) = Ini_Read_Int('sn1')
  sn(2) = Ini_Read_Int('sn2')
  simn = sn(2) - sn(1) + 1
  !//////////////////!

  file_stem = trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))//'_interp'//trim(RealToStr(interp_factor,3))//'_method'//trim(IntToStr(interp_method))//'_'

  if (want_pol) file_stem=trim(file_stem)//'pol_'
  file_stem = trim(file_stem)//trim(IntToStr(lens_method)) 
  cls_lensed_file  = trim(file_stem)//'.dat'
  
  call SetIdlePriority();

!Healpix 2.01 has up to 8192 
  if (w8name=='') then
    write (*,*) 'Warning: using unit weights as no w8dir specified'
    call HealpixInit(H,nside,lmax,.true.,w8dir='',method=mpi_division_method) 
  else
    call HealpixInit(H,nside,lmax,.true.,w8dir=w8name,method=mpi_division_method) 
  end if

  if (H%MpiID ==0) then !if we are main thread
    !All but main thread stay in HealpixInit
    call HealpixPower_nullify(P)
    call HealpixAlm_nullify(A)
    call HealpixMap_nullify(GradPhi)
    call HealpixMap_nullify(M)

    write(*,*) "Read Cls"
    !Reads in unlensed Cl text files as produced by CAMB
    if (TBEB) then
      f2 = Ini_Read_String('cls_phi')
      f3 = Ini_Read_String('cls_vio')
      call HealpixPower_ReadFromTextFile_EXT(P,f1,f2,f3,lmax,pol=want_pol,dolens=.true.)
    else
      call HealpixPower_ReadFromTextFile(P,f1,lmax,pol=want_pol,dolens=.true.)
    end if

    allocate(Cl(simn,9,lmax),Grad(0:lmax,0:lmax));  Cl=0d0;  Grad=0d0
    !!$OMP PARAllEl DO !PRIVATE(Cl)
    do n = sn(1), sn(2)
      ii = n - sn(1) + 1
      call MapSimulation(n,ii,Cl,Grad)
    end do
    !!$OMP END PARAllEl DO
    deallocate(Grad)

    !#### rlz averaged Cls ####!
    if(WantLCl) then
      allocate(mCl(10,lmax))
      mCl(1,:) = linspace(1,lmax)
      do i = 1, 9
        mCl(i+1,:) = sum(Cl(:,i,:),dim=1)/dble(simn)
      end do
      call savetxt('ucl.dat',mCl,ow=.true.)
      deallocate(mCl)
    end if
    !Free memory 
    deallocate(Cl)
    call HealpixAlm_Free(A)
   end if

  call Ini_Close

#ifdef MPIPIX
  call HealpixFree(H)
  call mpi_finalize(i)
#endif

#ifdef DEBUG
  write (*,*) 'End of program'
  pause
#endif

contains


subroutine MapSimulation(n,ii,Cl,Grad)
  implicit none
  !I/O
  integer, intent(in) :: n, ii
  real(dp), intent(out), optional :: Cl(:,:,:)
  complex(SPC), intent(inout) :: Grad(:,:)
  !internal
  character(len=80), dimension(1:60) :: hdr
  integer i 
  double precision, dimension(:,:), allocatable :: lCl
  real(SP), allocatable :: map(:,:)
  complex(SPC), allocatable :: alm(:,:,:)

  call write_minimal_header(hdr,'MAP',nside=nside,Ordering='RING')

  if(lensfix.and.n==1) then 
    write(*,*) "lensfix"
    call HealpixAlm_Sim(A,P,rand_seed,HasPhi=.true.,dopol = want_pol)
    Grad(2:lmax,0:lmax) = A%Phi(1,2:lmax,0:lmax)
  else if(lensfix.and.n>=2) then 
    write(*,*) "lensfix"
    call HealpixAlm_Sim(A,P,rand_seed,HasPhi=.true.,dopol = want_pol)
    A%Phi(1,2:lmax,0:lmax) = Grad(2:lmax,0:lmax)
  else
    if (.not.TBEB) call HealpixAlm_Sim(A,P,rand_seed,HasPhi=.true.,dopol = want_pol)
    if (TBEB)      call HealpixAlm_Sim_TEB(A,P,rand_seed)
  end if

  !* check cls
  if (present(Cl)) then
    call calccl(A%TEB(1,:,:),A%TEB(1,:,:),[1,lmax],Cl(ii,1,:))
    call calccl(A%TEB(2,:,:),A%TEB(2,:,:),[1,lmax],Cl(ii,2,:))
    call calccl(A%TEB(3,:,:),A%TEB(3,:,:),[1,lmax],Cl(ii,3,:))
    call calccl(A%TEB(1,:,:),A%TEB(2,:,:),[1,lmax],Cl(ii,4,:))
    call calccl(A%TEB(1,:,:),A%TEB(3,:,:),[1,lmax],Cl(ii,5,:))
    call calccl(A%TEB(2,:,:),A%TEB(3,:,:),[1,lmax],Cl(ii,6,:))
    call calccl(A%phi(1,:,:),A%Phi(1,:,:),[1,lmax],Cl(ii,7,:))
    call calccl(A%TEB(1,:,:),A%Phi(1,:,:),[1,lmax],Cl(ii,8,:))
    call calccl(A%TEB(2,:,:),A%Phi(1,:,:),[1,lmax],Cl(ii,9,:))
  end if

  if (WantGMAP) then 
    allocate(map(0:npix-1,1))
    write(*,*) "glm to gmap"
    call scalalm2map(H,A%lmax,A%Phi,map(:,1))
    call output_map(map,hdr,trim(out_file_root)//"Gmap_r"//str(n)//".fits")
    deallocate(map)
  end if

  if (Ini_Read_Logical('want_ucmb')) then
    allocate(map(0:npix-1,1))
    write(*,*) "ulm to ucmb"
    call scalalm2map(H,A%lmax,A%TEB,map(:,1))
    call output_map(map,hdr,trim(out_file_root)//'uTmap_r'//str(n)//'.fits')
    deallocate(map)
  end if

  write(*,*) "Alm gradient maps"
  call HealpixAlm2GradientMap(H,A, GradPhi,npix,'PHI')
  select case(lens_method)
  case (lens_exact)
    write(*,*) "Remapping (exact)"
    call HealpixExactLensedMap_GradPhi(H,A,GradPhi,M)
  case (lens_interp)
    write(*,*) "Remapping (interp)"
    call HealpixInterpLensedMap_GradPhi(H,A,GradPhi,M,interp_factor,interp_method)
  case default
    stop 'unknown lens_method'
  end select

  !#### compute Cl for each realization ####!
  if (WantLCl) then
    write(*,*) "Map to Alm"
    if (.not.want_pol)  call map2scalalm(H, lmax, M%TQU(:,1), A%TEB,-1.d0)
    if (want_pol)       call map2polalm(H, lmax, M%TQU, A%TEB,-1.d0)
    allocate(lCl(5,0:lmax))
    lcl(1,:) = linspace(1,lmax)
    do i = 1, 3
      call calccl(A%TEB(i,:,:),A%TEB(i,:,:),[1,lmax],lCl(i+1,:))
    end do
    call calccl(A%Phi(1,:,:),A%Phi(1,:,:),[1,lmax],lCl(5,:))
    call savetxt('lcl.dat',lcl,ow=.true.)
    deallocate(lCl)
  end if

  !#### Save map to fits files ####!
  write(*,*) "output lensed map"
  call output_map(real(M%TQU(0:npix-1,1:1)),hdr,trim(out_file_root)//'Tmap_r'//str(n)//'.fits')
  if (want_pol) then
    call output_map(real(M%TQU(0:npix-1,2:2)),hdr,trim(out_file_root)//"Qmap_r"//str(n)//".fits")
    call output_map(real(M%TQU(0:npix-1,3:3)),hdr,trim(out_file_root)//"Umap_r"//str(n)//".fits")
  end if
  write(*,*) "end lensed map"

end subroutine MapSimulation


subroutine HealpixPower_ReadFromTextFile_EXT(P,f1,f2,f3,lmax,pol,dolens) !namikawa
  use AMLutils
  Type(HealpixPower) P
  character(LEN=*), intent(IN) :: f1,f2,f3
  integer, intent(in) :: lmax
  logical,intent(in), optional :: dolens, pol
  integer i, l
  real(DP) :: T, E, B, TE, TB, EB, scal, phi, phiT, phiE

  call HealpixPower_Free(P)
  call HealpixPower_Init(P, lmax, pol = .true., dolens = .true.)

  open(unit=20,file=f1,status='old')
  open(unit=21,file=f2,status='old')
  open(unit=22,file=f3,status='old')
  P%Cl = 0
  P%PhiCl = 0
  TB=0d0; EB=0d0
  do i = 2, lmax
    read(21,*) l, T, E, B, TE, phi, phiT, phiE
    read(20,*) l, T, E, B, TE
    read(22,*) l, TB, EB
    if (l <= lmax) then
      scal = twopi/dble(l*(l+1))
      P%Cl(l,1) = T*scal
      P%Cl(l,2) = E*scal
      P%Cl(l,3) = B*scal
      P%Cl(l,4) = TE*scal
      P%Cl(l,5) = TB*scal
      P%Cl(l,6) = EB*scal
      P%PhiCl(l,1) = phi * scal/dble(l*(l+1))
      P%PhiCl(l,2) = phiT * scal/dble(l*(l+1))**0.5
      P%PhiCl(l,3) = phiE * scal/dble(l*(l+1))**0.5
    end if
  end do
  close(20)
  close(21)
  close(22)

end subroutine HealpixPower_ReadFromTextFile_EXT


end program SimLensCMB

