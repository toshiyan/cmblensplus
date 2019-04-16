! Simple program demonstrating how to generate a simulated lensed map
! AL, Feb 2004; Updated Oct 2007

! Modified by Toshiya Namikawa 

program SimLensCMB

  use HealpixObj
  use HealpixVis
  use Random
  use spinalm_tools
  use IniFile
  use AMLUtils

  !namikawa
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
  logical :: WantGMAP, WantLCl
  integer :: n, simn, el, i
  double precision, parameter :: Tcmb = 2.72d6, pi0 = 3.14d0
  double precision, dimension(:,:), allocatable :: mCl,vCl
  double precision, dimension(:,:,:), allocatable :: Cl
  !

  Ini_Fail_On_Not_Found = .true.
  call Ini_Open(GetParam(1), 3,err)
  if (err) stop 'No ini'
  nside  = Ini_Read_Int('nside')
  npix = 12*nside**2
  lmax   = Ini_Read_Int('lmax')  
  WantLCl = Ini_Read_Logical('WantLCl')
  f1 = Ini_Read_String('cls_file')
  f2 = Ini_Read_String('cls_phi')
  f3 = Ini_Read_String('cls_vio')
  out_file_root = Ini_Read_String('out_file_root')
  lens_method = Ini_Read_Int('lens_method')
  want_pol = Ini_Read_Logical('want_pol')
  rand_seed = Ini_Read_Int('rand_seed')
  interp_method = Ini_read_int('interp_method')
  Ini_Fail_On_Not_Found = .false.
  w8name = Ini_Read_String('w8dir')
  interp_factor=0
  if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
  simn = Ini_Read_Int('simn')
 
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
    write(*,*) "read cls"
    !Reads in unlensed Cl text files as produced by CAMB
    call HealpixPower_ReadFromTextFile(P,f1,lmax,pol=.true.,dolens=.true.)
    !call HealpixPower_ReadFromTextFile_EXT(P,f1,f2,f3,lmax,pol=.true.,dolens=.true.)

    allocate(Cl(simn,9,lmax));  Cl=0d0
    do n = 1, simn
      call MapSimulation(n,Cl)
    end do

    !#### averaged Cls ####!
    allocate(mCl(9,lmax),vCl(9,lmax))
    do el = 2, lmax
      do i = 1, 9
        mCl(i,el) = sum(Cl(:,i,el))/dble(simn)
      end do
    end do
    open(unit=20,file="mcl.dat",status="replace")
    do el = 2, lmax
      write(20,"(I6,9(1X,E14.7))") el, mCl(1:9,el)
    end do
    close(20)
    deallocate(mCl,vCl,Cl)
    call HealpixAlm_Free(A)
   end if

  call Ini_Close

contains


subroutine MapSimulation(n,Cl)
  implicit none
  !I/O
  integer, intent(in) :: n
  real(dp), intent(out), optional :: Cl(:,:,:)
  !internal
  character(len=80), dimension(1:60) :: hdr
  character(LEN=64) :: f, bname
  integer l
  real(8), dimension(:,:), allocatable :: cls

  call write_minimal_header(hdr,'MAP',nside=nside,Ordering='RING')

  write(bname,*) n
  call HealpixAlm_Sim(A,P,rand_seed,HasPhi=.true.,dopol = want_pol)

  !* check cls
  if (present(Cl)) then
    call calccl(A%TEB(1,:,:),A%TEB(1,:,:),[2,lmax],Cl(n,1,:))
    call calccl(A%TEB(2,:,:),A%TEB(2,:,:),[2,lmax],Cl(n,2,:))
    call calccl(A%TEB(3,:,:),A%TEB(3,:,:),[2,lmax],Cl(n,3,:))
    call calccl(A%TEB(1,:,:),A%TEB(2,:,:),[2,lmax],Cl(n,4,:))
    call calccl(A%TEB(1,:,:),A%TEB(3,:,:),[2,lmax],Cl(n,5,:))
    call calccl(A%TEB(2,:,:),A%TEB(3,:,:),[2,lmax],Cl(n,6,:))
    call calccl(A%phi(1,:,:),A%Phi(1,:,:),[2,lmax],Cl(n,7,:))
    call calccl(A%TEB(1,:,:),A%Phi(1,:,:),[2,lmax],Cl(n,8,:))
    call calccl(A%TEB(2,:,:),A%Phi(1,:,:),[2,lmax],Cl(n,9,:))
  end if

  write(*,*) "Alm gradient maps"
  call HealpixAlm2GradientMap(H,A, GradPhi,npix,'PHI')
  if (lens_method == lens_exact) then
    write(*,*) "Remapping (exact)"
    call HealpixExactLensedMap_GradPhi(H,A,GradPhi,M)
  else if (lens_method == lens_interp) then
    write(*,*) "Remapping (interp)"
    call HealpixInterpLensedMap_GradPhi(H,A,GradPhi,M,interp_factor,interp_method)
  else
    stop 'unknown lens_method'
  end if

  !#### compute Cl for each realization ####!
  write(*,*) "Map to Alm"
  call map2polalm(H, lmax, M%TQU, A%TEB,-1.d0)
  allocate(cls(4,0:lmax))
  call calccl(A%TEB(1,:,:),A%TEB(1,:,:),[2,lmax],cls(1,:))
  call calccl(A%TEB(2,:,:),A%TEB(2,:,:),[2,lmax],cls(2,:))
  call calccl(A%TEB(3,:,:),A%TEB(3,:,:),[2,lmax],cls(3,:))
  call calccl(A%Phi(1,:,:),A%Phi(2,:,:),[2,lmax],cls(4,:))
  open(unit=20,file="lensedCls.dat",status="replace")
  do l = 2, lmax
    write(20,"(I6,X,4(E14.7,X))") l, cls(1:4,l)
  end do 
  close(20)
  deallocate(cls)

  !#### Save map to fits files ####!
  write(*,*) "output observed map"
  call output_map(real(M%TQU(0:npix-1,1:1)),hdr,trim(out_file_root)//"Tmap_r"//trim(adjustl(bname))//".fits")
  call output_map(real(M%TQU(0:npix-1,2:2)),hdr,trim(out_file_root)//"Qmap_r"//trim(adjustl(bname))//".fits")
  call output_map(real(M%TQU(0:npix-1,3:3)),hdr,trim(out_file_root)//"Umap_r"//trim(adjustl(bname))//".fits")
  write(*,*) "end observed map"

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


subroutine CALCCL(alm1,alm2,el,Cl)
  implicit none
  !I/O
  integer, intent(in) :: el(2)
  double precision, intent(out) :: Cl(:)
  complex(spc), intent(in), dimension(0:el(2),0:el(2)) :: alm1, alm2
  !intenral
  integer :: l
  real(dp) :: tCl(el(2))

  tCl = 0d0
  do l = el(1), el(2)
    tCl(l) = ( REAL(alm1(l,0)*alm2(l,0)) + sum(alm1(l,1:l)*conjg(alm2(l,1:l))) + sum(alm2(l,1:l)*conjg(alm1(l,1:l))) )/(2.*l+1.)
  end do
  Cl = tCl

end subroutine CALCCL


end program SimLensCMB

