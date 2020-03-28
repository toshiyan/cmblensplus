!MPI Healpix routines, including generation of spin-s maps, 
!mapping gradients of scalars and the exact and approx weak lensed CMB 
!Antony Lewis 2004-2011, Based on Healpix 2

!Requires linking to Healpix libraries:
!See http://www.eso.org/science/healpix/

!Sign conventions follow Healpix/CMBFAST/CAMB
!Most easily used using HealpixObj.f90 wrapper routines
!Compile with -DMPIPIX -to use MPI

!Performance could be improved by changing the theta-CPU sharing
!However scaling is quite good up to 50 or so processors for high res transforms
!Temporary arrays use more memory than non-MPI routines
!For compatibility with Healpix input/output alm arrays are not packed (2*waste of mem)

!Jan 2005: improved/fixed polarization lens rotation factors. Minor fixes.
!Sept 2005: fixed bug in map2polalm
!Nov 2007: added bicubic interpolation, temp only, speedups
!Dec 2007: multiple map transforms, reduced memory requirements
!Jan 2008: further memory reductions for non-lensed; one MPI thread workarounds for scalar
!Oct 2010: corrected approximate handling of pole region interpolation (tiny area, virtually no effect)
!Nov 2010: fixes for bugs that only showed up in gfortran (thanks to Giancarlo de Gasperis)
!Apr 2011: Fixed wrap-around of phi during interp lensing

module MPIstuff
implicit none
double precision starttime
#ifdef MPIPIX
    include "mpif.h"
    integer ::  DebugMsgs =1
    integer MPIstatus(MPI_STATUS_SIZE), ierr
    integer SP_MPI,CSP_MPI 
#endif
contains
   
  subroutine MpiBarrier
  integer i
#ifdef MPIPIX
  call MPI_BARRIER(MPI_COMM_WORLD,i)
#endif
  end subroutine MpiBarrier

  subroutine GetMpiStat(MpiId, MpiSize)
   implicit none
   integer MpiId,MpiSize  
#ifdef MPIPIX  
   integer ierror

        call mpi_comm_rank(mpi_comm_world,MpiId,ierror)
        if (ierror/=MPI_SUCCESS) stop 'GetMpiDetail: MPI rank'
        call mpi_comm_size(mpi_comm_world,MpiSize,ierror)
        SP_MPI = MPI_REAL
        CSP_MPI=  MPI_COMPLEX
#else
  MpiId=0
  MpiSize=1   
#endif
  end subroutine GetMpiStat

 subroutine SyncInts(i,j,k)
  integer, intent(inout) :: i
  integer, intent(inout), optional :: j,k
#ifdef MPIPIX 

  integer params(3),sz

  params(1)=i
  sz=1
  if (present(j)) then
   params(2)=j
   sz=2
  end if
  if (present(k)) then
   params(3)=k
   sz=3
  end if

  call MPI_BCAST(params,sz,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  i= params(1)
  if (present(j)) then
   j=params(2)
  end if
  if (present(k)) then
   k=params(3)
  end if
    
#endif 
 end subroutine SyncInts

 subroutine SyncReals(i,j,k)
  real, intent(inout) :: i
  real, intent(inout), optional :: j,k
#ifdef MPIPIX 

  real params(3)
  integer sz

  params(1)=i
  sz=1
  if (present(j)) then
   params(2)=j
   sz=2
  end if
  if (present(k)) then
   params(3)=k
   sz=3
  end if

  call MPI_BCAST(params,sz,MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
  i= params(1)
  if (present(j)) then
   j=params(2)
  end if
  if (present(k)) then
   k=params(3)
  end if
    
#endif 
 end subroutine SyncReals

end module MPIstuff

module spinalm_tools
  use utilities, only: die_alloc
  use healpix_types
  use healpix_fft, only : real_fft
  IMPLICIT none

  Type HealpixInfo
     integer :: nside, lmax, Lastlmax
     logical pol
     REAL(KIND=DP), DIMENSION(:,:), Pointer :: w8ring_TQU  => NULL() 
     INTEGER(I8B), DIMENSION(:), pointer :: istart_south  => NULL() , istart_north => NULL()   
     COMPLEX(DPC),DIMENSION(:), pointer :: trig  => NULL() 
     REAL(DP), DIMENSION(:), Pointer :: recfac  => NULL() , Lambda_slm  => NULL() 
     integer MpiId, MPISize, MpiStat, last_nph
     integer(I4B), dimension(:), pointer :: ith_start => NULL() , ith_end => NULL() 
     integer, dimension(:), pointer :: North_Start => NULL() , North_Size => NULL() , &
        South_Start => NULL() , South_Size => NULL() 
  end type HealpixInfo

  Type HealpixMapArray
    REAL(SP), DIMENSION(:,:), pointer :: M  => NULL() 
  end Type HealpixMapArray

 type HealpixAllCl
  !All (a^i a^j) C_l
  !Index 0:lmax, i, j, where i,j are T E B
    real(SP), dimension(:,:,:), pointer :: Cl  => NULL() 
 end type HealpixAllCl 
  
 type HealpixCrossPowers
  !Array of cross-power spectra
  integer nmaps, lmax, npol 
  Type(HealpixAllCl), dimension(:,:), pointer :: Ps  => NULL() 
 end type HealpixCrossPowers

 type HealpixPackedScalAlms
   COMPLEX(SPC), dimension(:,:), pointer :: alms  => NULL() 
 end type HealpixPackedScalAlms

 type HealpixPackedAlms
   COMPLEX(SPC), dimension(:,:,:), pointer :: alms  => NULL() 
 end type HealpixPackedAlms
  

  Type LensGradients
     COMPLEX(SPC), DIMENSION(:), pointer :: grad_phiN => NULL() , grad_phiS => NULL() 
  end  Type LensGradients  

  integer, parameter :: interp_edge = 2
   !number of high-res pixels to go outside deflected region to get good interpolation
  
  integer, parameter :: EB_sign = -1
    !definition: for pol {}_2 a_{lm} = EB_sign*(E + iB)_{lm}
    !EB_sign = -1 corresponds to Healpix and CAMB/CMBFAST conventions

  logical :: mmax_approx = .true.
   !when true, uses that fact that don't need high m near the poles because the legendre 
   !functions are tiny for m >> l sin(theta)

  integer, parameter :: interp_basic=0, interp_cyl = 1

  integer, parameter :: division_equalrows=1, division_equalpix=2, division_balanced =3


  ! keep everything private unless stated otherwise
  private
  ! define large and small numbers used to renormalise the recursion on the Legendre Polynomials
  real(KIND=DP), private, PARAMETER :: FL_LARGE = 1.0e30_dp
  real(KIND=DP), private, PARAMETER :: FL_SMALL = 1.0e-30_dp
  real(KIND=DP), private :: OVFLOW, UNFLOW, ScaleFactors(-10:10)  
 
  ! make (front end) routines public
  public :: spinalm2map, alm2GradientMap, map2spinalm,scalalm2map, mmax_approx, HealpixInfo, &
            HealpixInit,HealpixFree, map2scalalm, a_ix, scalalm2LensedMap, &
            alm2Lensedmap, map2polalm, polalm2map, alm2LensedQuadContrib, EB_sign, &
            alm2LensedmapInterp, scalalm2LensedmapInterp,scalalm2LensedmapInterpCyl, &
            alm2LensedmapInterpCyl, interp_basic, interp_cyl , GeteTime, &
            division_equalrows, division_equalpix, division_balanced, HealpixMapArray, &
            HealpixCrossPowers, HealpixAllCl, maparray2scalcrosspowers,maparray2crosspowers, &
            HealpixCrossPowers_Free, healpix_wakeMPI, healpix_sleepMPI, scalalm2bispectrum
contains

 function GeteTime()
      use MPIStuff
      double precision GeteTime
#ifndef MPIPIX
      real etime

      call cpu_time(etime)    
      GeteTime = etime
#else
      GeteTime = MPI_WTime()
#endif
 end function GeteTime


  function a_ix(lmax, l, m) result(index)
    integer, intent(in) :: lmax, l, m
    integer :: index
    index = (m*(2*lmax-m+1))/2 + l + 1
  end function a_ix

   subroutine HealpixInit(H, nside, lmax, HasPol, w8dir, method)
    USE fitstools, ONLY : getsize_fits, input_map
    use MPIStuff
    use healpix_types
    Type (HealpixInfo) :: H
    Integer, optional, intent(in) :: method
    logical, intent(in), optional :: HasPol
    character(LEN=*), optional, intent(in) :: w8dir
    logical use_weights
    integer, intent(in) :: nside, lmax
    !real(dp) logOVFLOW
    integer npol, n_rings
    character(LEN=120) :: sstr, filename
    REAL(SP), DIMENSION(:,:), allocatable :: w8
    integer nph, i, delta, st
    ! Changed for new division between threads
    Real(sp) :: mean_pix !Mean number of pixels in each section of northern hemisphere
    real(sp) :: pix_w, mean_weight,time_weights(2*nside)
    Integer ::  pixels, row
    Integer :: division = 1 
      !Determines whether to give each section
      ! equal numbers of rows (1), or equal numbers of pixels (2)
      ! (2) is much faster for 'exact' lensing
#ifdef MPIPIX    
    integer status, ierror
#endif  
    CHARACTER(LEN=*), PARAMETER :: code = 'HealpixInit'
 

#ifndef MPIPIX
        call HealpixFree(H)
!If MPI must call healpixFree manually
#endif    
       nullify(H%recfac,H%Lambda_slm)

       call HealpixInitTrig(H,nside,lmax)

       npol = 1
       if (present(HasPol)) then
        npol = 3
       end if
       H%pol = HasPol
  

    allocate(H%w8ring_TQU(1:2*nside,1:max(1,npol)))
    
    use_weights = present(w8dir) 
    if (use_weights) then
     use_weights = w8dir/='' 
    end if
    if (use_weights) then
         allocate(w8(1:2*nside,1:max(1,npol)))
 
         write (sstr,"(I5.5)") nside
         filename= trim(w8dir)//"weight_ring_n"//trim(sstr)//".fits"

         n_rings = 2 * nside
     
         if (getsize_fits(filename) /= n_rings) then
            write (*,*) 'HealpixInit:wrong file'//trim(filename)
            stop
         endif
     
         if (HasPol) then
            call input_map(filename, w8, n_rings, 3, fmissval=0.0_sp)
         else
            call input_map(filename, w8, n_rings, 1, fmissval=0.0_sp)
         endif

         H%w8ring_TQU =  1 + w8
         deallocate(w8)
     else

       H%w8ring_TQU=1
  
     endif

!Get factors for making well behaved Ylm
    OVFLOW=exp(log(FL_LARGE))
    UNFLOW=exp(log(FL_SMALL))
   ! logOVFLOW=log(FL_LARGE)
    ScaleFactors=0
    do i=-10,10
     ScaleFactors(i) = FL_LARGE**i !exp(i*logOVFLOW)
    end do

!Mpi properties
    H%MpiId = 0; H%MpiSize = 1
    H%MpiStat = 0


#ifdef MPIPIX
        if (SP==KIND(1.d0)) then
         SP_MPI = MPI_DOUBLE_PRECISION
         CSP_MPI = MPI_DOUBLE_COMPLEX
        else if (SP == KIND(1.)) then
         SP_MPI = MPI_REAL
         CSP_MPI=  MPI_COMPLEX
        else
         stop 'Unknown SP KIND for MPI'
        end if
        
        call mpi_comm_size(mpi_comm_world,H%MpiSize,ierror)
        call mpi_comm_rank(mpi_comm_world,H%MpiId,ierror)
        if (ierror/=MPI_SUCCESS) stop 'HealpixInit: MPI rank'
#endif

!Sectioning of the sphere between threads
!Following things to bear in mind:
! * range of l needed smaller near poles
! * healpix has 4*nside pix per ring for i>nside, but linear with i for i<=nside
! * - this means some rings have inefficient FFT at i<nside where i is still large
! * map2alm and alm2map are naively ~ proportional only to number of rings if FFT efficient
! * Lensing interpolation time is roughly proportional to number of pixels

        If (present(method)) division = method
#ifdef MPIPIX      
        if (DebugMsgs >1 .and. H%MpiId==0) print *,'mpi_division = ', division
        if (H%MpiSIze==1) division = division_equalrows 
#else
        division = division_equalrows 
#endif

!       If (division == division_balanced) Then
!        st = (2*nside)/(3*H%MpiSize)      !Put more in poles for balanced
!        delta = (2*nside - st)/H%MpiSize
!        st = 1 + st + mod(2*nside-st,H%MpiSize)
!       else 
        delta = (2*nside)/H%MpiSize
        st = 1 + mod(2*nside,H%MpiSize)
!       end if
       
       allocate(H%ith_start(0:H%MpiSIze-1), H%ith_end(0:H%MpiSIze-1), H%North_Start(0:H%MpiSIze-1), &
         H%North_Size(0:H%MpiSIze-1), H%South_Start(0:H%MpiSIze-1), H%South_Size(0:H%MpiSIze-1))
       H%ith_start = 1
       H%ith_end = 2*nside

       if ( division == division_balanced) then
              
        do i=1, nside*2
         !for healpix transform timing grows approximately linearly to nside, then drops and stays ~constant 
          if (i< nside) then
              time_weights(i) = 0.7 +  i*24./nside     !
          else
            time_weights(i) = 16 !+ real(2*(i-nside))/nside
          end if 
!          time_weights(i) = 8 + nside*(sin(i*pi/(4*nside))**0.8 + 0.2)  + max(0,(i-nside)/10)
!!          if (i> nside/2 .and. i< nside) then
!           time_weights(i) = time_weights(i)*1.2        
!          end if     

        end do
         mean_weight = sum(time_weights)/H%MpiSize
         pix_w = 0
        
       else if  (division == division_equalpix)  then 

!Giving less to poles usually a good idea
!        if (H%MpiSize<3) then
!         first_pix = (3*nside*(6*nside+2)/H%MpiSize)/4
!        else
!         first_pix = (nside*(6*nside+2)/H%MpiSize)/2
!        end if
!        mean_pix = (nside*(6*nside+2) -  first_pix)/(H%MpiSize-1)

        pixels = 0
        mean_pix = nside*(6*nside+2)/H%MpiSize
       end if


        do i= 0, H%MpiSize -1 
    
           If (division == division_equalpix) Then


              ! New version SJS 15/12/2004 for equal pixels per thread
              ! mean_pix is the average number of pixels given to each thread

              ! New method - divide into ~equal numbers of pixels
              ! Should be significantly faster if using 'exact' lensing method
              ! very marginaly faster if using interpolation method 
              ! ideally need a third method which gives less to the poles
              ! if doing interpolation
              if (i == 0) then
                 H%ith_start(i) = 1
              else
                 H%ith_start(i) = H%ith_end(i-1) + 1
              end if

              row = H%ith_start(i)-1
              do while (pixels .LT. (i+1.0)*mean_pix)
                 row = row + 1
                 nph = 4*nside
                 if (row .LT. nside) nph = 4*row
                 pixels = pixels + nph
              end do
              H%ith_end(i) = row
              If (i == (H%MpiSize-1)) H%ith_end(i) = 2*nside

           Else If (division == division_equalrows) Then
              !divide into equal numbers of rows

            if (i == 0) then
             !Do more, but poles are faster anyway
              H%ith_start(i) = 1
              H%ith_end(i) = st + delta-1 
            else
              H%ith_start(i) = st + i*delta
              H%ith_end(i) =  H%ith_start(i) +  delta -1
            end if
            
           Else If (division == division_balanced) Then
            !New method Oct 07
            
              if (i == 0) then
                 H%ith_start(i) = 1
              else
                 H%ith_start(i) = H%ith_end(i-1) + 1
              end if

              row = H%ith_start(i)-1
              do while (pix_w < (i+1)*mean_weight .and. row < 2*nside)
                 row = row + 1
                 pix_w = pix_w + time_weights(row)
              end do
              If (i == (H%MpiSize-1)) row= 2*nside
              H%ith_end(i) = row
#ifdef MPIPIX
          if (DebugMsgs > 1 .and. H%MpiId==0) write(*,*) i, 'row end = ',row      
#endif            
           Else
              Stop 'HealpixInit : Unknown method'
           End If


            if (H%ith_end(i)< nside) then  
                  nph = 4*H%ith_end(i)
               else                  
                  nph = 4*nside
            endif
       
            H%North_start(i) = H%istart_north(H%ith_start(i)-1)
            H%North_Size(i) = H%istart_north(H%ith_end(i)-1) + nph &
                              -H%North_start(i)
        
            if (H%ith_start(i) < nside) then  
                  nph = 4*H%ith_start(i)
               else                   
                  nph = 4*nside
            endif
            if (H%ith_end(i) == nside*2) then
              H%South_start(i) = H%istart_south(H%ith_end(i)-1)
            else
             H%South_start(i) = H%istart_south(H%ith_end(i))
            end if
            H%South_Size(i) = H%istart_south(H%ith_start(i)) + nph &
                              - H%South_start(i)

        end do
#ifdef MPIPIX
        if (H%MpiId>0) call MessageLoop(H)
#endif
 
  end subroutine HealpixInit

   subroutine HealpixInitTrig(H, nside, lmax, not_healpix)
    use MPIStuff
    use healpix_types
    logical, intent(in), optional :: not_healpix
    logical not_heal
    Type (HealpixInfo) :: H
    integer, intent(in) :: lmax, nside
    integer ith, status, nph
    CHARACTER(LEN=*), PARAMETER :: code = 'HealpixTrig'


       nullify(H%trig)
       H%last_nph = -1
       H%lmax = lmax
   !    H%nalms_max = ((lmax+1)*(lmax+2))/2
       H%Lastlmax = 0
       H%nside = nside
      
      not_heal = .false.
      if (present(not_healpix)) not_heal = not_healpix
      
      if (not_heal) then
       nullify(H%istart_north)
       nullify(H%istart_south)
      else
      
        ALLOCATE(H%istart_north(0:2*nside),stat = status)
        if (status /= 0) call die_alloc(code,'istart_north')

        ALLOCATE(H%istart_south(0:2*nside),stat = status)
        if (status /= 0) call die_alloc(code,'istart_south')

        H%istart_north(0)=0
        H%istart_south(0)=12*int(nside,I8B)**2
        do ith=1,2*nside
           if (ith.lt.nside) then  ! polar cap (north)
              nph = 4*ith
           else                   ! tropical band (north) + equator
              nph = 4*nside
           endif
           H%istart_north(ith)=H%istart_north(ith-1)+nph
           H%istart_south(ith)=H%istart_south(ith-1)-nph
        enddo
        
      end if  

   end subroutine HealpixInitTrig

   subroutine HealpixInfo_GetTrig(H, nph)
     Type (HealpixInfo) :: H
     integer, intent(in) :: nph    
     integer status, m
     real(dp) phi0    

     if (H%last_nph /= nph) then
 
       deallocate(H%trig,stat = status)
       ALLOCATE(H%trig(0:max(2*H%nside,H%lmax)),stat = status) 
 
       H%trig=1
       phi0=PI/DBLE((nph/4)*4)
       do m=0,max(2*H%nside,H%lmax)
          H%trig(m)= CMPLX( DCOS(m*phi0), DSIN(m*phi0), kind=DP)
       enddo
       H%last_nph = nph
      end if

   end  subroutine HealpixInfo_GetTrig

  
   function NearestFastFFTnum(i)
    !returns next number of form 2^n 3^m for low m
    integer, intent(in) :: i
    integer NearestFastFFTnum
    integer j, vals(71)
    
    vals = (/128, 144, 192, 256, 288, 384, 512, 576, 768, 1024, 1152, 1536, 2048, 2304, 3072, & 
     4096, 4608, 6144, 8192, 9216, 12288, 16384, 18432, 24576, 32768, 36864, 49152, 65536, 73728, &
      98304, 131072, 147456, 196608, 262144, 294912, 393216, 524288, 589824, 786432, 1048576, &
      1179648, 1572864, 2097152, 2359296, 3145728, 4194304, 4718592, 6291456, 8388608, &
      9437184, 12582912, 16777216, 18874368, 25165824, 33554432, 37748736, 50331648, 67108864, &
       75497472, 100663296, 134217728, 150994944, 201326592, 268435456, 301989888, 402653184, 452984832, &
       536870912, 603979776, 805306368, 905969664/)

    do j=1,71
     if (i*0.9<=vals(j)) then
      NearestFastFFTnum = vals(j)
      return
     end if
    end do
   
    stop 'NearestFastFFTnum: number too large'
   
   end function NearestFastFFTnum
  
   function ScaleFactor(i)
    integer, intent(in) :: i
    real(dp) :: ScaleFactor
    
     if (i>-10) then
        ScaleFactor = ScaleFactors(i)
     else
        ScaleFactor = 0
     end if

   end function ScaleFactor

   subroutine HealpixFree(H)
    Type (HealpixInfo) :: H
    integer status
#ifdef MPIPIX
    if (H%MpiId == 0) call SendMessages(H, 'EXIT')
#endif    
    deallocate(H%w8ring_TQU, stat = status)
    deallocate(H%istart_north, H%istart_south, stat = status)
    deallocate(H%trig, stat = status)
    deallocate(H%recfac, stat = status)
    deallocate(H%ith_start, H%ith_end,H%North_Start, H%North_Size, & 
           H%South_Start, H%South_Size, stat = status)
    nullify(H%w8ring_TQU)
   end subroutine HealpixFree


   subroutine HealpixInitRecfac(H,nlmax)
     Type (HealpixInfo) :: H
     INTEGER(I4B), intent(in):: nlmax
     integer(I8B) :: m, l  
     integer status, a_ix
     integer l2, m2

     if (H%MpiId > 0 .and. associated(H%recfac) .and. nlmax == H%Lastlmax) return   
     call HealpixFreeRecfac(H)  
     H%Lastlmax = nlmax
     deallocate(H%recfac,stat= status)
     ALLOCATE(H%recfac(((nlmax+1)*(nlmax+2))/2),stat = status)    
     if (status /= 0) call die_alloc('HealpixInitRecfac','recfac')

     a_ix = 0
     do m = 0, nlmax
      m2 = m**2
      do l = m, nlmax
        a_ix = a_ix + 1
        l2 = (l+1)**2
        H%recfac(a_ix) = SQRT( real(4 * l2 - 1,dp) / real(l2-m2,dp) )
      end do
     end do 
 
   end subroutine HealpixInitRecfac

   
          
  subroutine HealpixFreeRecfac(H)
      Type (HealpixInfo) :: H
         integer status

         if (H%MpiId > 0) return   !cache it as have loads of memory

         deallocate(H%recfac,stat= status)


  end subroutine HealpixFreeRecfac


   function get_mmax(nlmax,sth)
        integer, intent(in) :: nlmax
        real(dp), intent(in) :: sth
        integer get_mmax
        
        if (mmax_approx) then
            get_mmax = min(nlmax,max(40,nint(1.25*nlmax*sth)))
        else
            get_mmax = nlmax
        end if

   end function get_mmax


  function l_min_ylm(m, sth) result(lmin)
  !From heapix 2, roughly consistent with choice of get_mmax above
  !================================================================
    ! returns minimal order l at which to keep Ylm
    ! |Ylm| < eps * Y00 ==>
    ! m_cut(theta, l) = theta * l * e / 2 + | ln(eps)| + ln(l)/2
    ! if eps = 1.e-15 and l < 1.e4
    ! m_cut(theta, l) = theta * l * 1.35 + 40
    ! the choice of 1.35 (or larger) 
    ! also insures that the equatorial rings will have all their Ylm's computed
    ! default parameters are HPX_MXL0 = 40 and HPX_MXL1 = 1.35_DP
    !======================================================
    ! parameters of short-cut: defined in module header
    ! dummy variables
    integer(I4B)             :: lmin
    integer(I4B), intent(IN) :: m
    real(DP),     intent(IN) :: sth
    integer, parameter :: HPX_MXL0 = 40
    real(dp), parameter :: HPX_MXL1 = 1.35_dp
    
    lmin = m ! default
    if (mmax_approx) lmin = max(lmin, int((m - HPX_MXL0)/(HPX_MXL1 * sth)))

    return
  end function l_min_ylm
 
  subroutine spinring_synthesis(H,nlmax,datain,nph,dataout,kphi0,mmax_ring)
 !Don't fully follow the signs here, but the answer is correct
 !Note no point using FFTW etc as FFT is a negligible fraction of computation cost
    !=======================================================================
    !     RING_SYNTHESIS
    !       called by alm2map
    !       calls     real_fft
    !
    !     dataout(j) = sum_m datain(m) * exp(i*m*phi(j)) 
    !     with phi(j) = j*2pi/nph + kphi0*pi/nph and kphi0 =0 or 1
    !
    !     as the set of frequencies {m} is larger than nph, 
    !     we wrap frequencies within {0..nph-1}
    !     ie  m = k*nph + m' with m' in {0..nph-1}
    !     then
    !     noting bw(m') = exp(i*m'*phi0) 
    !                   * sum_k (datain(k*nph+m') exp(i*k*pi*kphi0))
    !        with bw(nph-m') = CONJ(bw(m')) (if datain(-m) = CONJ(datain(m)))
    !     dataout(j) = sum_m' [ bw(m') exp (i*j*m'*2pi/nph) ]
    !                = Fourier Transform of bw
    !        is real
    !
    !         NB nph is not necessarily a power of 2
    !
    !=======================================================================

    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax
    INTEGER(I4B), INTENT(IN) :: mmax_ring
    INTEGER(I4B), INTENT(IN) :: nph, kphi0

    COMPLEX(DPC), DIMENSION(0:nlmax), INTENT(IN) :: datain
    REAL(SP),     DIMENSION(0:nph-1), INTENT(OUT)     :: dataout
    REAL(DP),     DIMENSION(0:nph-1)     :: data
    INTEGER(I4B) :: iw,ksign,m,k,kshift
    COMPLEX(DPC), DIMENSION(0:nph-1) :: bw
    COMPLEX(DPC) :: dat
#ifdef MPIPIX    
    integer status
#endif  
   
    !=======================================================================

    call HealpixInfo_GetTrig(H, nph)

    nsmax = H%nside
    ksign = + 1

    kshift = (-1)**kphi0  ! either 1 or -1
    bw(0:nph-1) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

    !     all frequencies [-m,m] are wrapped in [0,nph-1]
    bw(0)=datain(0)
    do m  = 1, mmax_ring                        ! in -nlmax, nlmax
       iw = MODULO(m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + datain(m)*(kshift**k)  ! complex number
       iw = MODULO(-m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (-m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + CONJG(datain(m))*(kshift**k)  ! complex number
    enddo
    !     kshift**k = 1       for even turn numbers
    !               = 1 or -1 for odd  turn numbers : results from the shift in space

    !     applies the shift in position <-> phase factor in Fourier space
    data(0)=REAL(bw(0))

    !Data is in packed storage
    do iw = 1, nph/2  -1
       m = ksign*(iw)
       if(kphi0==1) then
          dat =bw(iw) * H%trig(m)
       else
          dat =bw(iw)
       endif
       data(iw*2-1 ) = REAL(dat)
       data(iw*2) = AIMAG(dat)

    enddo
!    nph is always even for Healpix
    iw=nph/2
   m = ksign*(iw)
   if(kphi0==1) then
       dat =bw(iw) * H%trig(m)
    else
      dat =bw(iw)
    endif
    data(iw*2-1) = REAL(dat)
 
    call real_fft (data, backward=.true.)
    !     ^^^^^^^^^^^^
    dataout=REAL(data(0:nph-1))

    RETURN
  END subroutine spinring_synthesis
  
  subroutine alm2GradientMap(H, inlmax, alm, map_QU)
 !Get the map of the gradient of alm (know pure E, so quicker than general routine)
 !internally use EB_sign=1 convention, though result is independent
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H
    INTEGER(I4B), INTENT(IN) :: inlmax 
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:)  :: alm
    COMPLEX(SPC), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map_QU
    COMPLEX(SPC), DIMENSION(:), pointer :: map2N,map2S
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0, nlmax

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, lambda_w_1, lambda_x_1, a_w 
    REAL(DP) :: a_w_m
    COMPLEX(DPC) :: factor_1, factor_2, factor_1_1, factor_2_1
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2GRADIENTMAP'
    COMPLEX(DPC) ::  b_north_Q(0:H%lmax), b_north_U(0:H%lmax)
    COMPLEX(DPC) ::  b_south_Q(0:H%lmax), b_south_U(0:H%lmax)
    INTEGER(I4B) :: status,par_lm, a_ix
    REAL(DP) , DIMENSION(:), ALLOCATABLE :: cth_l
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:),   ALLOCATABLE :: ringR, ringI
    integer mmax_ring, nalms, lmin
#ifdef MPIPIX    
    double precision Initime
#endif      
    !=======================================================================

     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call SyncInts(nlmax)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms),stat = status )
     if (status /= 0) call die_alloc(code,'alm2')
     if (H%MpiId==0) call Alm2PackAlm(alm,alm2,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code//' Got alm ',H%MpiId, GeteTime() - StartTime
     allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1),stat = status) 
     if (status /= 0) call die_alloc(code,'map2')   
     allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1),stat = status) 
     if (status /= 0) call die_alloc(code,'map2')   
#else
     map2N => map_QU
     map2S => map_QU
#endif

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(cth_l(nlmax))    

    ALLOCATE(ringR(0:4*nsmax-1),ringI(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    call HealpixInitRecfac(H,nlmax)
    call GetLamfact(lam_fact, nlmax)
   
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )
    !     --------------------------------------------

    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)   ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       do l=1, nlmax
        cth_l(l) = cth*real(l,dp)
       end do
       
       mmax_ring = get_mmax(nlmax,sth)
            
       b_north_Q(0:nlmax) = 0
       b_north_U(0:nlmax) = 0
       b_south_Q(0:nlmax) = 0
       b_south_U(0:nlmax) = 0
       
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1
       a_ix = 0
       a_w = -1._dp / sth

       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = -1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW

          lam_lm = corfac     !  actual lambda_mm      
          
          a_ix = a_ix + 1
          if (m >=1) then
!normal_l cancels with gradient, sign from gradient 
              lambda_x = - lam_lm * fm / sth 
              lambda_w = -lambda_x * cth

              b_n_Q =  lambda_w * alm2(a_ix)
              b_s_Q =  par_lm * b_n_Q

              b_n_U = (0,-1)* lambda_x * alm2(a_ix)
              b_s_U = -par_lm * b_n_U

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          
          lmin  = l_min_ylm(m, sth)
!         
!          do l = m+1, nlmax
!             par_lm = - par_lm  ! = (-1)^(l+m+s)
!             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
!             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
!            
!             a_ix = a_ix + 1
!
!             if (l >= lmin) then
!              lambda_x = a_w * fm  * lam_lm
!              lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - real(l,dp)*cth*lam_lm) 
!
!              factor_1 =  lambda_w * alm2(a_ix)
!              b_n_Q = b_n_Q +          factor_1 
!              b_s_Q = b_s_Q + par_lm * factor_1 ! X has a diff. parity
!
!              factor_2 =   (0,1) *  lambda_x * alm2(a_ix) 
!              b_n_U = b_n_U - factor_2
!              b_s_U = b_s_U + par_lm * factor_2
!             end if
!            
!             lam_0 = lam_1 / a_rec
!             lam_1 = lam_2
!             a_rec = H%recfac(a_ix)
!             lam_2 = (cth * lam_1 - lam_0) * a_rec
!
!             if (abs(lam_2)  >  OVFLOW) then
!                lam_0=lam_0/OVFLOW
!                lam_1=lam_1/OVFLOW
!                lam_2 = (cth * lam_1 - lam_0) * a_rec
!                scalel=scalel+1
!                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!             elseif (abs(lam_2)  <  UNFLOW) then
!                lam_0=lam_0*OVFLOW
!                lam_1=lam_1*OVFLOW
!                lam_2 = (cth * lam_1 - lam_0) * a_rec 
!                scalel=scalel-1
!                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!             endif
!
!          enddo
          a_w_m = a_w*fm 
          do l = m+1, nlmax-1, 2
          !This is semi-optimized version where we do two at once
          !par_lm starts off positive  (negative on entry to loop) 
          !doesn't gain much here       

             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
           
             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_ix = a_ix + 1
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec
       
             if (l >= lmin) then

              lambda_w_1 = a_w * (lam_fact(a_ix)*lam_lm1m - cth_l(l)*lam_lm) 
              lambda_x_1 = a_w_m * lam_lm

              lam_lm1m=lam_lm  ! actual lambda_l-1,m 
              lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed

              lambda_w = a_w * (lam_fact(a_ix+1)*lam_lm1m - cth_l(l+1)*lam_lm) 
              lambda_x = a_w_m  * lam_lm

              factor_1_1 =  lambda_w_1 * alm2(a_ix)
              factor_1 =  lambda_w * alm2(a_ix+1)
       
              b_n_Q = b_n_Q +  factor_1 + factor_1_1
              b_s_Q = b_s_Q -  factor_1 + factor_1_1

              factor_2_1 = lambda_x_1*cmplx(-aimag(alm2(a_ix)),real(alm2(a_ix)))
              factor_2 =   lambda_x*cmplx(-aimag(alm2(a_ix+1)),real(alm2(a_ix+1))) 

              b_n_U = b_n_U - factor_2 - factor_2_1
              b_s_U = b_s_U - factor_2 + factor_2_1
             else
              lam_lm1m=lam_lm  ! actual lambda_l-1,m 
              lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             end if
             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_ix = a_ix + 1
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2+lam_1)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2+lam_1)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo
          if (mod(nlmax - m,2)==1) then
          !do last one

             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed

             a_ix = a_ix + 1
            
             lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - cth_l(nlmax)*lam_lm) 
             lambda_x = a_w * fm  * lam_lm

             factor_1 =  lambda_w * alm2(a_ix)
             b_n_Q = b_n_Q +  factor_1 
             b_s_Q = b_s_Q +  factor_1  

             factor_2 =   (0,1) *  lambda_x * alm2(a_ix) 
             b_n_U = b_n_U - factor_2
             b_s_U = b_s_U + factor_2
     
!           print *,b_n_Q,b_n_U,b_s_Q,b_s_U 
            
          end if

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U

       enddo

       call spinring_synthesis(H,nlmax, b_north_Q, nph, ringR, kphi0,mmax_ring)
       call spinring_synthesis(H,nlmax, b_north_U, nph, ringI, kphi0,mmax_ring)
       map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))
       if (ith  <  2*nsmax) then
          call spinring_synthesis(H,nlmax, b_south_Q, nph, ringR, kphi0,mmax_ring)
          call spinring_synthesis(H,nlmax, b_south_U, nph, ringI, kphi0,mmax_ring)
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))
       endif

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecfac(H)
    DEALLOCATE(lam_fact, cth_l)
    DEALLOCATE(ringR,ringI)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiId
    StartTime = Getetime()
    call MPI_GATHERV(map2N(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),CSP_MPI, &
       map_QU,H%North_Size,H%North_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2S(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),CSP_MPI, &
       map_QU,H%South_Size,H%South_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *,code //' Wait at gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2N,map2S) 
#endif

  end subroutine alm2GradientMap

  subroutine GetLamfact(lam_fact, nlmax)
  real(dp) lam_fact(*), fm2
  integer, intent(in) :: nlmax
  integer a_ix,l,m

    a_ix = 0
    do m = 0, nlmax
      fm2 = real(m,dp) **2
      a_ix = a_ix + 1
      do l = m+1, nlmax
        a_ix = a_ix + 1
        lam_fact(a_ix) = SQRT( (2 * l + 1) / real(2*l - 1,dp) * (l**2-fm2))
      end do
   end do 

  end subroutine GetLamFact


  subroutine spinalm2map(H,inlmax, alm_EB, map_QU, inspin)
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B), INTENT(IN) :: inlmax, inspin
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_EB
    COMPLEX(SPC), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map_QU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: EB
    COMPLEX(SPC), DIMENSION(:), pointer :: map2
    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) ::  factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U

    CHARACTER(LEN=*), PARAMETER :: code = 'SPINALM2MAP'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax, spin

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:),   ALLOCATABLE :: ringR, ringI
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer a_ix, nalms
#ifdef MPIPIX
    double precision Initime
#endif
    !=======================================================================

    !     --- allocates space for arrays ---


     nsmax = H%nside
     nlmax = inlmax
     spin = inspin

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
     call SyncInts(nlmax,spin)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(EB(2,nalms))
     if (H%MpiId==0) call EB2PackEB(alm_EB,EB,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(EB,SIze(EB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)    
     if (status /= 0) call die_alloc(code,'map2')
#else
     map2 => map_QU 
#endif

    if (spin<1 .or. spin>3) stop 'Only spin 1 to 3 supported'

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')

    ALLOCATE(ringR(0:4*nsmax-1),ringI(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')


    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   if (spin ==2 ) lam_fact = lam_fact * 2 !HealPix polarization def
   

    normal_l = 0.0_dp
    do l = spin, nlmax
       fl = DBLE(l)
       if (spin==1) then 
        normal_l(l) = EB_sign*SQRT( 1 / ((fl+1.0_dp)*fl) ) 
       else if (spin==2) then
        normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
       else if (spin==3) then
        normal_l(l) = EB_sign*SQRT( 1 / ((fl+3.0_dp)*(fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)*(fl-2.0_dp)) ) 
       end if 
    enddo

    !     ----- set the whole map to zero ------
!    map_QU = 0.0
    !     --------------------------------------
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !     --------------------------------------------

    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       mmax_ring = get_mmax(nlmax,sth) 
       
       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = (-1)**spin  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac     !  actual lambda_mm      
          a_ix = a_ix + 1
          !l=m special case
          if (m >=spin) then

              if (spin==1) then
                !swapped
                lambda_x = normal_l(m) * lam_lm * fm / sth 
                lambda_w = -lambda_x * cth
              else if (spin==2) then
                lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
                lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              else if (spin==3) then
                lambda_x = normal_l(m) * lam_lm / sth * fm*(fm-1)*(fm-2) 
                lambda_w = lambda_x * cth * ( 1 - 4*one_on_s2)
                lambda_x = lambda_x * (4*one_on_s2 - 3)             
              end if 
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
              b_n_Q =  lambda_w * EB(1,a_ix) + zi_lam_x * EB(2,a_ix)
              b_s_Q =  par_lm*(lambda_w * EB(1,a_ix) - zi_lam_x * EB(2,a_ix))

              ! - alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
              b_n_U = lambda_w * EB(2,a_ix) - zi_lam_x * EB(1,a_ix)
              b_s_U = par_lm*(lambda_w * EB(2,a_ix) + zi_lam_x * EB(1,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1
             if (l>=spin .and. corfac /= 0) then

                 if (spin==1) then
                    a_w = normal_l(l) / sth
                    lambda_x = a_w * fm  * lam_lm
                    lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 
                 else if (spin==2) then
                     a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                     b_w =  c_on_s2 * lam_fact(a_ix)
                     a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                     lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                     lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 else if (spin==3) then
                     a_w = normal_l(l) /sth
                     b_w = (l-1)*(l-2)
                     lambda_x = a_w*fm*( (one_on_s2*(4*fm2-(12*l-8)) - 3 * b_w)*lam_lm + &
                          12*lam_fact(a_ix) * c_on_s2  * lam_lm1m)  
                     lambda_w = a_w*( (l*b_w - one_on_s2*(8*l+fm2*(4*l-12))) * cth * lam_lm  &
                        - lam_fact(a_ix)*( fl2 + fl +6 - (8+4*fm2)*one_on_s2) * lam_lm1m)

                 end if
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * EB(1,a_ix)
                 factor_2 =  zi_lam_x * EB(2,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * EB(2,a_ix) 
                 factor_2 =   zi_lam_x * EB(1,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U

       enddo

       call spinring_synthesis(H,nlmax, b_north_Q, nph, ringR, kphi0,mmax_ring)
       call spinring_synthesis(H,nlmax, b_north_U, nph, ringI, kphi0,mmax_ring)
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))
  
       if (ith  <  2*nsmax) then
          call spinring_synthesis(H,nlmax, b_south_Q, nph, ringR, kphi0,mmax_ring)
          call spinring_synthesis(H,nlmax, b_south_U, nph, ringI, kphi0,mmax_ring)
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))

       endif

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)
    DEALLOCATE(ringR,ringI)
    deallocate(EB)

#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiId,' Time: ', GeteTime() - iniTime
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),CSP_MPI, &
       map_QU,H%North_Size,H%North_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),CSP_MPI, &
       map_QU,H%South_Size,H%South_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
!    call MPI_REDUCE(map2,map,size(map),SP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine spinalm2map

  subroutine map2spinalm(H, inlmax, map_QU,alm_EB, spinin, cos_theta_cut)
    use MPIStuff
    Type (HealpixInfo) :: H

    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    COMPLEX(SPC), INTENT(OUT),  DIMENSION(:,:,:) :: alm_EB
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: map_QU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: EB
    COMPLEX(SPC), DIMENSION(:), pointer :: map2

    integer, intent(in) :: spinin
    REAL(DP),     INTENT(IN) :: cos_theta_cut

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax, spin       ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x

    CHARACTER(LEN=*), PARAMETER :: code = 'MAP2SPINALM'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_nQ, phas_nU
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_sQ, phas_sU
    INTEGER(I4B) mmax_ring, status, par_lm, a_ix
    integer nalms

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    LOGICAL   :: keep_it
#ifdef MPIPIX
    double precision Initime
#endif
    !=======================================================================

    nsmax = H%nside
    nlmax = inlmax
    spin = spinin
     
#ifdef MPIPIX
     if (cos_theta_cut/=-1) stop 'cos_theta_cut /= -1'
     if (H%MpiId==0) then 
      if(DebugMsgs>0) print *,code //': Sending to farm '
      call SendMessages(H,code)
      map2 => map_QU
    else
       allocate(map2(0:12*H%nside**2-1),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
    end if

    StartTime = getetime()    
    iniTime = StartTime
    call SyncInts(nlmax,spin)

    call MPI_SCATTERV(map_QU,H%North_Size, H%North_Start, &
       CSP_MPI, map2(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_SCATTERV(map_QU,H%South_Size, H%South_Start, &
       CSP_MPI, map2(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *,code //' Scattered ',H%MpiId, GeteTime() - StartTime
#else
    map2 => map_QU
#endif

    nalms = ((nlmax+1)*(nlmax+2))/2   

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(phas_nQ(0:nlmax),&
         &   phas_nU(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_sQ(0:nlmax),&
         &   phas_sU(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    if (spin<1 .or. spin>3) stop 'Only spin 1 to 3 supported'

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamfact(lam_fact, nlmax)
   if (spin==2) lam_fact = lam_fact*2

   allocate(EB(2,nalms),stat = status)
   if (status /= 0) call die_alloc(code,'EB')
   EB = 0    
       
   omega_pix = EB_sign * pi / (3 * nsmax * real(nsmax,dp))

    normal_l = 0.0_dp
    do l = spin, nlmax
       fl = DBLE(l)
       if (spin==1) then 
        normal_l(l) = omega_pix * SQRT( 1 / ((fl+1.0_dp)*fl) ) 
       else if (spin==2) then
        normal_l(l) = omega_pix * SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
       else if (spin==3) then
        normal_l(l) = omega_pix * SQRT( 1 / ((fl+3.0_dp)*(fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)*(fl-2.0_dp)) ) 
       end if 
    enddo

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    
    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)

       phas_nQ=0; phas_sQ=0;phas_nU=0;phas_sU=0

       if (ith  <=  nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2

       mmax_ring = get_mmax(nlmax,sth) 

       keep_it = (ABS(cth) > cos_theta_cut) ! part of the sky out of the symmetric cut
       if (keep_it) then
          ring(0:nph-1) = real(map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_nQ, kphi0, mmax_ring)
          ring(0:nph-1) = aimag(map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_nU, kphi0, mmax_ring)
       endif

       if (ith  <  2*nsmax .and. keep_it) then
          ring(0:nph-1) = real(map2(H%istart_south(ith):H%istart_south(ith)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_sQ, kphi0, mmax_ring)
          ring(0:nph-1) = aimag(map2(H%istart_south(ith):H%istart_south(ith)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_sU, kphi0, mmax_ring)
       endif

       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       if (keep_it) then ! avoid un-necessary calculations (EH, 09-2001)
          lam_mm = sq4pi_inv
          scalem=1
          a_ix = 0
          do m = 0, mmax_ring
             fm  = DBLE(m)
             f2m = 2.0_dp * fm
             fm2 = fm * fm
             fm_on_s2 = fm * one_on_s2

             !           ---------- l = m ----------
             par_lm = (-1)**spin   ! = (-1)^(l+m+s)
             if (m  >=  1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm) < UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif

             a_ix = a_ix+1

             corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
             lam_lm = corfac
        
            if (m >=spin) then
            
              if (spin==1) then
                lambda_x = normal_l(m) * lam_lm * fm / sth 
                lambda_w = -lambda_x * cth
              else if (spin==2) then
                lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
                lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              else if (spin==3) then
                lambda_x = normal_l(m) * lam_lm / sth * fm*(fm-1)*(fm-2) 
                lambda_w = lambda_x * cth * ( 1 - 4*one_on_s2)
                lambda_x = lambda_x * (4*one_on_s2 - 3)             
              end if 
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
              
                 EB(1,a_ix) = EB(1,a_ix) &
                  &                 + lambda_w * (phas_nQ(m) + par_lm*phas_sQ(m)) &
                  &                 + zi_lam_x * (phas_nU(m) - par_lm* phas_sU(m))

                 EB(2,a_ix) = EB(2,a_ix) &
                  &                 + lambda_w * (phas_nU(m) + par_lm*phas_sU(m)) &
                  &                 - zi_lam_x * (phas_nQ(m) - par_lm*phas_sQ(m))
            
            end if

             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec
             do l = m+1, nlmax
                par_lm = - par_lm  ! = (-1)^(l+m)
                lam_lm1m=lam_lm ! actual lambda_l-1,m (useful for polarisation)
                lam_lm   = lam_2*corfac ! actual lambda_lm (OVFLOW factors removed)
                fl  = DBLE(l)
                fl2 = fl * fl

                a_ix = a_ix + 1
             if (l>=spin .and. corfac /= 0) then
                 !Corfac=0 guarantees lam(l-1) is also v close to zero
                  
                 if (spin==1) then
                    a_w = normal_l(l) / sth
                    lambda_x = a_w * fm  * lam_lm
                    lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 
                 else if (spin==2) then
                     a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                     b_w =  c_on_s2 * lam_fact(a_ix)
                     a_x =  2*(l-1)* cth * lam_lm
                     lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                     lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 else if (spin==3) then
                     a_w = normal_l(l) /sth
                     b_w = (l-1)*(l-2)
                     lambda_x = a_w*fm*( (one_on_s2*(4*fm2-(12*l-8)) - 3 * b_w)*lam_lm + &
                          12*lam_fact(a_ix) * c_on_s2  * lam_lm1m)  
                     lambda_w = a_w*( (l*b_w - one_on_s2*(8*l+fm2*(4*l-12))) * cth * lam_lm  &
                        - lam_fact(a_ix)*( fl2 + fl +6 - (8+4*fm2)*one_on_s2) * lam_lm1m)

                 end if

                zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 EB(1,a_ix) = EB(1,a_ix) &
                     &          + lambda_w * (phas_nQ(m) + par_lm*phas_sQ(m)) &
                     &          + zi_lam_x * (phas_nU(m) - par_lm*phas_sU(m))
                 EB(2,a_ix) = EB(2,a_ix) &
                     &         +  lambda_w * (phas_nU(m) + par_lm*phas_sU(m)) &
                     &         - zi_lam_x  * (phas_nQ(m) - par_lm*phas_sQ(m))

              end if ! l allowed by spin or zero

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_2)  >  OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                elseif (abs(lam_2)  <  UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                endif

             enddo ! loop on l
          enddo ! loop on m
       endif ! test on cut sky
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(phas_nQ,phas_nU)
    DEALLOCATE(phas_sQ,phas_sU)
    DEALLOCATE(ring)
#ifdef MPIPIX
    if (H%MpiId>0) deallocate(map2)
    StartTime = Getetime()
    if (H%MpiId==0) then
     call MPI_REDUCE(MPI_IN_PLACE,EB,size(EB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    else
     call MPI_REDUCE(EB,MPI_IN_PLACE,size(EB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    end if
    if (DebugMsgs>1) print *,code//' done reduce ', H%MpiId, GeteTime() -StartTime
    if (H%MpiId == 0) call PackEB2EB(EB,alm_EB, nlmax)
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code //' Time: ',GeteTime() - IniTime
#else
    call PackEB2EB(EB,alm_EB, nlmax)
#endif
    deallocate(EB)

  END subroutine map2spinalm

  subroutine spinring_analysis(H, nlmax, datain,nph,dataout,kphi0,mmax_ring)
    !=======================================================================
    !     ring_analysis
    !       called by map2alm
    !       calls     real_fft
    !
    !     integrates (data * phi-dependence-of-Ylm) over phi
    !     --> function of m can be computed by FFT
    !     with  0<= m <= npoints/2 (: Nyquist)
    !     because the data is real the negative m are the conjugate of the 
    !     positive ones
    !=======================================================================
    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax
    INTEGER(I4B), INTENT(IN) :: mmax_ring
    INTEGER(I4B), INTENT(IN) :: nph, kphi0

    REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
    COMPLEX(DPC), DIMENSION(0:nlmax), INTENT(OUT) :: dataout
    INTEGER(I4B) :: i,m,im_max,ksign
    REAL(DP), DIMENSION(0:nph-1) :: data
#ifdef MPIPIX    
    integer status
#endif  
    !-----------------------------------------------------------------------

    call HealpixInfo_GetTrig(H,nph)

    nsmax = H%nside
 
    ksign = - 1
    data=0.
    data(0:nph-1)=datain(0:nph-1)

    call real_fft(data, backward=.false.)

    im_max = MIN(nph/2,mmax_ring)
    dataout(0)=CMPLX(data(0),0.0_dp,kind=DP)

    do i = 1, im_max*2-3, 2
       dataout((i+1)/2) = CMPLX( data(i), data(i+1),kind= DP) 
    enddo

    if(im_max==nph/2) then
       dataout(im_max)= CMPLX( data(nph-1),0.0_dp,kind=DP)
    else
       dataout(im_max)= CMPLX( data(2*im_max-1),data(2*im_max),kind=DP)
    endif

    if(im_max==mmax_ring) goto 1000

    do i =  im_max+1,min(nph,mmax_ring)
       dataout(i) = conjg(dataout(2*im_max-i) )
    end do

    if(min(nph,mmax_ring)==mmax_ring) goto 1000

    do i =  2*im_max+1,mmax_ring
       dataout(i) = dataout(mod(i,2*im_max)) 
    end do

1000 continue

    if(kphi0==1)then
       do i =  0,mmax_ring
          m = ksign*i
          dataout(i)=dataout(i)* CONJG(H%trig(-m))
       enddo
    end if


  END subroutine spinring_analysis


    !=======================================================================
  subroutine scalalm2map(H, inlmax, alm, map)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature field
    !     map(theta,phi) = sum_l_m a_lm Y_lm(theta,phi)
    !                    = sum_m {e^(i*m*phi) sum_l a_lm*lambda_lm(theta)}
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the sum over m is done by FFT
    !
    !              -------------------------------
    !          precomputes the Lambda_lm recurrence factor 
    !      and is therefore ~50% faster than previous versions
    !              -------------------------------
    !
    !=======================================================================
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
    REAL(SP),     INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map

    REAL(SP),     DIMENSION(:), pointer :: map2N, map2S
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax          ! alm related
    INTEGER(I4B) :: nph, kphi0                         ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2
    REAL(DP) :: f2m, corfac
    COMPLEX(DPC) :: b_n, b_s, factor, factor2

    CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2MAP'
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring !,  par_lm
    integer nalms, a_ix

    REAL(SP), DIMENSION(0:4*H%nside-1) :: ring
#ifdef MPIPIX    
    double precision Initime
    integer status
#endif      
!=======================================================================
  
     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call SyncInts(nlmax)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms))
     if (H%MpiId==0) call Alm2PackAlm(alm,alm2,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
     if (H%MpiId==0) then
      map2N => map
      map2S => map
     else
      allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1),stat = status) 
      if (status /= 0) call die_alloc(code,'map2')   
      allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1),stat = status) 
      if (status /= 0) call die_alloc(code,'map2')
     end if   
#else
     map2S => map 
     map2N => map 
#endif

    call HealpixInitRecfac(H,nlmax)
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith =H%ith_start(H%MpiId), H%ith_end(H%MpiId)  
       !        cos(theta) in the pixelisation scheme

       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       mmax_ring = get_mmax(nlmax,sth) 

       lam_mm = sq4pi_inv ! lambda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          f2m = 2.0_dp * m

          !           ---------- l = m ----------
!          par_lm = 1  ! = (-1)^(l+m)
          if (m >= 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
  
          lam_lm = corfac
          a_ix = a_ix + 1
          b_n = lam_lm * alm2(a_ix)
          b_s = b_n

          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp 
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec

          do l = m+1, nlmax-1, 2
             
            a_ix = a_ix+1

            lam_0 = lam_1 / a_rec
            lam_1 = lam_2
    
            a_rec = H%recfac(a_ix)
            lam_2 = (cth * lam_1 - lam_0) * a_rec
            
            factor = (lam_1*corfac) * alm2(a_ix)
            factor2 = (lam_2*corfac) * alm2(a_ix+1)

            b_n = b_n + factor + factor2
            b_s = b_s - factor + factor2
            
            lam_0 = lam_1 / a_rec
            lam_1 = lam_2
            a_ix = a_ix+1
            a_rec = H%recfac(a_ix)
            lam_2 = (cth * lam_1 - lam_0) * a_rec
           
             if (abs(lam_1+lam_2) > OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_1+lam_2) < UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif
          enddo
          if (mod(nlmax-m,2)==1) then
                 a_ix = a_ix+1
                 b_n = b_n + corfac*lam_2*alm2(a_ix)
                 b_s = b_s - corfac*lam_2*alm2(a_ix)
          end if

!          do l = m+1, nlmax
!             par_lm = - par_lm  ! = (-1)^(l+m)
!
!             lam_lm = lam_2*corfac ! Remove OVFLOW-factors 
!             a_ix = a_ix + 1
!             factor = lam_lm * alm2(a_ix)
!             b_n = b_n +          factor
!             b_s = b_s + par_lm * factor
!
!             lam_0 = lam_1 / a_rec
!             lam_1 = lam_2
!             a_rec = H%recfac(a_ix)
!             lam_2 = (cth * lam_1 - lam_0) * a_rec
!
!             if (abs(lam_2) > OVFLOW) then
!                lam_0=lam_0/OVFLOW
!                lam_1=lam_1/OVFLOW
!                lam_2 = (cth * lam_1 - lam_0) * a_rec
!                scalel=scalel+1
!                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!             elseif (abs(lam_2) .lt. UNFLOW) then
!                lam_0=lam_0*OVFLOW
!                lam_1=lam_1*OVFLOW
!                lam_2 = (cth * lam_1 - lam_0) * a_rec 
!                scalel=scalel-1
!                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!             endif
!          enddo

          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H,nlmax,b_north,nph,ring,kphi0,mmax_ring)   ! north hemisph. + equator
       map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = ring(0:nph-1)
       
       if (ith < 2*nsmax) then
          call spinring_synthesis(H,nlmax,b_south,nph,ring,kphi0,mmax_ring) ! south hemisph. w/o equat
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1) = ring(0:nph-1)
       endif

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecFac(H)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code //' Gather ',H%MpiId
    
    StartTime = Getetime()
    if (H%MpiSize>1) then
    if (H%MpiID==0) then
     call MPI_GATHERV(MPI_IN_PLACE,H%North_Size(H%MpiId),SP_MPI, &
       map,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
     call MPI_GATHERV(MPI_IN_PLACE,H%South_Size(H%MpiId),SP_MPI, &
       map,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr) 
    else
     call MPI_GATHERV(map2N(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),SP_MPI, &
       map,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
     call MPI_GATHERV(map2S(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),SP_MPI, &
       map,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr) 
     if (H%MpiId > 0) deallocate(map2N,map2S)
    end if
    end if
    if (DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code //' Time :', GeteTime()-IniTime
    

#endif
  end subroutine scalalm2map


    subroutine map2scalalm(H,inlmax, map, alm, cos_theta_cut)
    !=======================================================================
    !     computes the a(l,m) from a map for the HEALPIX pixelisation
    !      for the Temperature field
    !     a(l,m) = int T(theta,phi) Y_lm(theta,phi)^* dtheta*dphi
    !            = int dtheta lambda_lm(theta) 
    !                  * int dphi T(theta,phi) e^(-i*m*phi)
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the integral over phi is done by FFT
    !
    !     cos_theta_cut (>0) is the cosine of the 
    !     symmetric cut applied to the sky
    !     if it is <= 0 no cut is applied
    !
    !     NB : these a(l,m) have to be multiplied by the pixel solid angle
    !      to give the correct coefficients
    !
    !             -------------------------------
    !         precomputes the Lambda_lm recurrence factor
    !      and is therefore ~50% faster than previous versions
    !     the multiplication by omega_pix is done in the routine
    !             -------------------------------
    !
    !=======================================================================
    use MPIStuff
    Type (HealpixInfo) :: H
    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    REAL(SP),     INTENT(IN),  DIMENSION(0:12*H%nside**2-1), target :: map
    COMPLEX(SPC), INTENT(OUT), DIMENSION(:,:,:) :: alm
    COMPLEX(SPC),   DIMENSION(:), allocatable :: alm2
    REAL(SP),     DIMENSION(:), pointer :: map2N,map2S
    
    REAL(DP),     INTENT(IN) :: cos_theta_cut

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax, a_ix       ! alm related
    INTEGER(I4B) :: nph, kphi0   ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_0, lam_1, lam_2
    REAL(DP) ::  f2m, corfac

    CHARACTER(LEN=*), PARAMETER :: code = 'MAP2SCALALM'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_n
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_s
    INTEGER(I4B) :: mmax_ring, status, nalms
    COMPLEX(DPC) phas_p, phas_m

    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    LOGICAL   :: keep_it
    integer lmin !, par_lm
#ifdef MPIPIX
    double precision Initime
#endif
    !=======================================================================


    nsmax = H%nside
    nlmax = inlmax
     
#ifdef MPIPIX
     if (H%MpiId==0) then 
      if(DebugMsgs>1) print *,code //': Sending to farm '
      IniTime = GeteTime()
      call SendMessages(H,code)
      map2N => map
      map2S => map
    else
       allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
       allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
    end if

    StartTime = getetime()    
    call SyncInts(nlmax)

     !Inplace doesn't seem to work here?? Also bug in openmpi 1.2.2 for MpiSize=1
    if (H%MpiSize>1) then
     call MPI_SCATTERV(map,H%North_Size, H%North_Start, &
       SP_MPI, map2N(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
     call MPI_SCATTERV(map,H%South_Size, H%South_Start, &
       SP_MPI, map2S(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end if
    if(DebugMsgs>1) print *,code //' Scattered ',H%MpiId, GeteTime() - StartTime
#else
    map2N => map
    map2S => map
#endif


    !     --- allocates space for arrays ---

    ALLOCATE(phas_n(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_s(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

    call HealpixInitRecfac(H,nlmax)

    nalms = ((nlmax+1)*(nlmax+2))/2   
    allocate(alm2(nalms), stat=status)
    if (status /= 0) call die_alloc(code,'alm2')

    alm2 = 0

    !     -------------------------------------------

    omega_pix = pi / (3.0_dp * nsmax * nsmax)

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)

       phas_n(0:nlmax) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)   ! North    m >= 0
       phas_s(0:nlmax) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)   ! South    m >= 0

       if (ith .le. nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif

       mmax_ring = get_mmax(nlmax,sth) 

       keep_it = (ABS(cth).gt.cos_theta_cut) ! part of the sky out of the symmetric cut


       if (keep_it) then
          ring(0:nph-1) = map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_n, kphi0, mmax_ring)
  
       if (ith .lt. 2*nsmax) then
          ring(0:nph-1) = map2S(H%istart_south(ith):H%istart_south(ith)+nph-1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_s, kphi0,mmax_ring)
       endif

       endif


       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       if (keep_it) then ! avoid un-necessary calculations (EH, 09-2001)

          lam_mm = sq4pi_inv * omega_pix ! lambda_00 * norm
          scalem=1 
          a_ix = 0
          do m = 0, mmax_ring
             f2m = 2.0_dp * m
         
             !           ---------- l = m ----------
             if (m .ge. 1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm).lt.UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif
             corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
             
             phas_p=phas_n(m) + phas_s(m)
             phas_m=phas_n(m) - phas_s(m)
 
             a_ix = a_ix + 1
             alm2(a_ix) = alm2(a_ix) + corfac * phas_p
             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec

             lmin = l_min_ylm(m, sth)
             do l = m+1, nlmax-1, 2
   
                a_ix = a_ix+1

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
      
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec
      
                if (l >= lmin) then
                 alm2(a_ix) = alm2(a_ix) + (lam_1*corfac) * phas_m
                 alm2(a_ix+1) = alm2(a_ix+1) + (lam_2*corfac) * phas_p
                end if

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_ix = a_ix+1
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_1+lam_2) .gt. OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                elseif (abs(lam_1+lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                endif
             enddo ! loop on l
             if (mod(nlmax-m,2)==1) then
                a_ix = a_ix+1
                alm2(a_ix) = alm2(a_ix) + lam_2*corfac * phas_m
             end if

!             par_lm=1
!             do l = m+1, nlmax
!                par_lm = - par_lm  ! = (-1)^(l+m)
!                
!   
!                lam_lm = lam_2*corfac ! Remove OVFLOW-factors
!                a_ix = a_ix+1
!                alm2(a_ix) = alm2(a_ix) + lam_lm * (phas_n(m) + par_lm*phas_s(m))
!
!                lam_0 = lam_1 / a_rec
!                lam_1 = lam_2
!                a_rec = H%recfac(a_ix)
!                lam_2 = (cth * lam_1 - lam_0) * a_rec
!
!                if (abs(lam_2) .gt. OVFLOW) then
!                   lam_0=lam_0/OVFLOW
!                   lam_1=lam_1/OVFLOW
!                   lam_2 = (cth * lam_1 - lam_0) * a_rec
!                   scalel=scalel+1
!                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!                elseif (abs(lam_2) .lt. UNFLOW) then
!                   lam_0=lam_0*OVFLOW
!                   lam_1=lam_1*OVFLOW
!                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
!                   scalel=scalel-1
!                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!                endif
!             enddo ! loop on l

!             if (m < nlmax) then
!             
!             lmin = l_min_ylm(m, sth)
!             lstart =m+1
!             if (ith==2*nsmax) then
!              !On equator Y_lm are zero if l+m is odd, need to make sure we get the rescaling
!              !so make sure phase of the pairing is correct. Note in this case phas_m=phas_p
!              lstart = lstart+1
!              lam_0 = lam_1 / a_rec
!              lam_1 = lam_2
!              a_ix = a_ix+1
!              a_rec = H%recfac(a_ix)
!              lam_2 = (cth * lam_1 - lam_0) * a_rec
!             end if
!             do l = lstart, nlmax-1, 2
!   
!                a_ix = a_ix+1
!
!                lam_0 = lam_1 / a_rec
!                lam_1 = lam_2
!      
!                a_rec = H%recfac(a_ix)
!                lam_2 = (cth * lam_1 - lam_0) * a_rec
!      
!                if (l >= lmin) then
!                 alm2(a_ix) = alm2(a_ix) + (lam_1*corfac) * phas_m
!                 alm2(a_ix+1) = alm2(a_ix+1) + (lam_2*corfac) * phas_p
!                end if
!
!                lam_0 = lam_1 / a_rec
!                lam_1 = lam_2
!                a_ix = a_ix+1
!                a_rec = H%recfac(a_ix)
!                lam_2 = (cth * lam_1 - lam_0) * a_rec
!
!                if (abs(lam_2) .gt. OVFLOW) then
!                   lam_0=lam_0/OVFLOW
!                   lam_1=lam_1/OVFLOW
!                   lam_2 = (cth * lam_1 - lam_0) * a_rec
!                   scalel=scalel+1
!                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!                elseif (abs(lam_2) .lt. UNFLOW) then
!                   lam_0=lam_0*OVFLOW
!                   lam_1=lam_1*OVFLOW
!                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
!                   scalel=scalel-1
!                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
!                endif
!             enddo ! loop on l
!             if (mod(nlmax-lstart,2)==0) then
!                a_ix = a_ix+1
!                alm2(a_ix) = alm2(a_ix) + lam_2*corfac * phas_m
!             end if
!             end if

          enddo ! loop on m
       endif ! test on cut sky
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(phas_n)
    DEALLOCATE(phas_s)
    DEALLOCATE(ring)
    call HealpixFreeRecFac(H)
    if (H%MpiId>0) deallocate(map2N,map2S)
#ifdef MPIPIX
    StartTime = Getetime()
    if (H%MpiId==0) then
     call MPI_REDUCE(MPI_IN_PLACE,alm2,size(alm2),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    else
     call MPI_REDUCE(alm2,MPI_IN_PLACE,size(alm2),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    end if
    if (DebugMsgs>1) print *,code //': Time at reduce ', H%MpiId, GeteTime() -StartTime
    if (H%MpiId == 0) then
      call PackAlm2Alm(alm2, alm, nlmax)
      print *,code // ' Time: ', GeteTime() - IniTime
    end if
#else
    call PackAlm2Alm(alm2, alm, nlmax)
#endif
    deallocate(alm2)

  END subroutine map2scalalm


  
  subroutine HealpixCrossPowers_Free(P)
   Type (HealpixCrossPowers) :: P
   integer  i, j
   
     do i=1, P%nmaps
      do j=1, i
       deallocate(P%Ps(i,j)%Cl)
       if (i/=j .and. P%npol>1) deallocate(P%Ps(j,i)%Cl)
      end do
     end do 
     deallocate(P%Ps)

  end subroutine HealpixCrossPowers_Free

  subroutine maparray2scalcrosspowers(H,inlmax, maps, P, nmaps, free)
    !=======================================================================
    ! Compute power array 1/(2l+1)\sum_m <alm(1)alm(2)^*> for all maps
    !=======================================================================
    use MPIStuff
    Type (HealpixInfo) :: H
    INTEGER(I4B), INTENT(IN) :: inlmax
    Type (HealpixMapArray), target :: maps(:)
    Type (HealpixCrossPowers) :: P
    integer, intent(in) :: nmaps
    logical, intent(in), optional :: free
    logical dofree
    
    Type(HealpixPackedScalAlms) :: A
    CHARACTER(LEN=*), PARAMETER :: code = 'MAPARRAY2SCALCROSSPOWERS'
    integer l,i,j,m, ix
    real(DP), dimension(:,:), allocatable :: Cl
#ifdef MPIPIX    
    double precision Initime
#endif      
    if (present(free)) then
     dofree=free
    else
     dofree = .false.
    end if  

    call maparray2packedscalalms(H,inlmax, maps, A, nmaps, dofree)

    if (H%MpiId==0) then
#ifdef MPIPIX
     if(DebugMsgs>1) print *,code //': Getting C_l '
     IniTime = GeteTime()
#endif
     P%nmaps = nmaps
     P%lmax = inlmax
     P%npol = 1
     allocate(Cl(nmaps,nmaps))
     allocate(P%Ps(nmaps,nmaps))
     do i=1, nmaps
      do j=1, i   
!       deallocate(P%Ps(i,j)%Cl,stat = status)
       allocate(P%Ps(i,j)%Cl(0:inlmax,1,1))
       if (i/=j) P%Ps(j,i)%Cl=>P%Ps(i,j)%Cl
      end do

     end do  
     
     do l=0, inlmax
          ix = a_ix(inlmax, l, 0)
          do i=1, nmaps
           do j=1, i   
             Cl(j,i)= REAL(A%alms(i,ix))*REAL(A%alms(j,ix))
            end do
           end do 
 
          do m=1, l
           ix = a_ix(inlmax, l, m)
              do i=1, nmaps
               do j=1, i   
                Cl(j,i) = Cl(j,i) + 2*REAL(A%alms(i,ix)*CONJG(A%alms(j,ix)))
               end do
             end do  
         end do !m

        do i=1, nmaps
         do j=1, i   
           P%Ps(i,j)%Cl(l,1,1) =Cl(j,i)/(2*l+1)
         end do
        end do 
     end do !l
     deallocate(Cl)

     deallocate(A%alms)
#ifdef MPIPIX
     if (DebugMsgs>0) print *,code // ' Time: ', GeteTime() - iniTime
#endif
    end if
    
  end subroutine maparray2scalcrosspowers   

  subroutine maparray2crosspowers(H,inlmax, maps, P, nmaps, free)
    !=======================================================================
    ! Compute power array 1/(2l+1)\sum_m <alm(1)alm(2)^*> for all maps
    !=======================================================================
    use MPIStuff
    Type (HealpixInfo) :: H
    INTEGER(I4B), INTENT(IN) :: inlmax
    Type (HealpixMapArray), target :: maps(:)
    Type (HealpixCrossPowers) :: P
    integer, intent(in) :: nmaps
    Type(HealpixPackedAlms) :: A
    Type(HealpixPackedScalAlms) :: AT
    CHARACTER(LEN=*), PARAMETER :: code = 'MAPARRAY2CROSSPOWERS'
    integer l,i,j,m, ix, polx,poly
    real(DP), dimension(:,:,:,:), allocatable :: Cl
    logical, intent(in), optional :: free
    logical dofree
    integer nT,nP
    Type (HealpixMapArray) :: MapsT(nmaps), MapsP(nmaps)
    integer numMaps(nmaps), indices(nmaps)

#ifdef MPIPIX    
    double precision Initime
#endif      
    if (present(free)) then
     dofree=free
    else
     dofree = .false.
    end if
    nT=0  
    do i=1,nmaps
     if (size(maps(i)%M,2)==1) then
      nT=nT+1
      numMaps(i)=1
      indices(i)=nT
      if (dofree) then
        allocate(mapsT(nT)%M(0:size(Maps(i)%M,1)-1,1))
        mapsT(nT)%M = Maps(i)%M
        deallocate(Maps(i)%M)
      else
        MapsT(nT)%M => maps(i)%M
      end if
     end if 
    end do
    
    if (nT==0) then !no maps with only temperature
     call maparray2packedpolalms(H,inlmax, maps, A, nmaps, dofree)
     numMaps=3
     do i=1,nmaps
      indices(i)=i
     end do
    else
        nP=0  
        do i=1,nmaps
         if (size(maps(i)%M,2)==3) then
          nP=nP+1
          numMaps(i)=3
          indices(i)=nP
          if (dofree) then
            allocate(mapsT(nP)%M(0:size(Maps(i)%M,1)-1,3))
            mapsP(nP)%M = Maps(i)%M
            deallocate(Maps(i)%M)
          else
            MapsP(nP)%M => maps(i)%M
          end if
         end if 
        end do
      call maparray2packedscalalms(H,inlmax, mapsT, AT, nT, dofree)
      call maparray2packedpolalms(H,inlmax, mapsP, A, nP, dofree)
    
    end if
    if (H%MpiId==0) then
#ifdef MPIPIX
     if(DebugMsgs>1) print *,code //': Getting C_l '
     IniTime = GeteTime()
#endif
     P%nmaps = nmaps
     P%lmax = inlmax
     P%npol = 3
     allocate(Cl(nmaps,nmaps,3,3))
     allocate(P%Ps(nmaps,nmaps))
     do i=1, nmaps
      do j=1, nmaps
      !We are duplicating work when polx==poly, never mind   
       allocate(P%Ps(i,j)%Cl(0:inlmax,numMaps(i),numMaps(j)))
       P%Ps(i,j)%Cl=0
      end do
     end do  
     
     do l=0, inlmax
          ix = a_ix(inlmax, l, 0)
          do polx=1,3
           do poly=1,polx
            do i=1, nmaps
             if (polx > numMaps(i)) cycle
             do j=1, nmaps   
             
              if (numMaps(j)==3 .and. numMaps(i)==3) then
                Cl(j,i,poly,polx)= REAL(A%alms(indices(i),polx,ix))*REAL(A%alms(indices(j),poly,ix))
              else 
               if (poly > numMaps(j)) cycle
               if (numMaps(j)==1 .and. numMaps(i)==3) then
                 Cl(j,i,poly,polx)= REAL(A%alms(indices(i),polx,ix))*REAL(AT%alms(indices(j),ix))
               else  if (numMaps(j)==3 .and. numMaps(i)==1) then
                 Cl(j,i,poly,polx)= REAL(AT%alms(indices(i),ix))*REAL(A%alms(indices(j),poly,ix))
               else
                 Cl(j,i,poly,polx)= REAL(AT%alms(indices(i),ix))*REAL(AT%alms(indices(j),ix))
               end if
              
               end if
             end do
            end do 
           end do
          end do 
         
          do m=1, l
           ix = a_ix(inlmax, l, m)
          
           do polx=1,3
            do poly=1,polx
              do i=1, nmaps
              if (polx > numMaps(i)) cycle
               do j=1, nmaps   
                  if (numMaps(j)==3 .and. numMaps(i)==3) then
                    Cl(j,i,poly,polx) = Cl(j,i,poly,polx) + 2*REAL(A%alms(indices(i),polx,ix)*CONJG(A%alms(indices(j),poly,ix)))
                  else 
                   if (poly > numMaps(j)) cycle
                   if (numMaps(j)==1 .and. numMaps(i)==3) then
                    Cl(j,i,poly,polx) = Cl(j,i,poly,polx) + 2*REAL(A%alms(indices(i),polx,ix)*CONJG(AT%alms(indices(j),ix)))
                   else  if (numMaps(j)==3 .and. numMaps(i)==1) then
                    Cl(j,i,poly,polx) = Cl(j,i,poly,polx) + 2*REAL(AT%alms(indices(i),ix)*CONJG(A%alms(indices(j),poly,ix)))
                   else
                    Cl(j,i,poly,polx) = Cl(j,i,poly,polx) + 2*REAL(AT%alms(indices(i),ix)*CONJG(AT%alms(indices(j),ix)))
                   end if
                  
                  end if

               end do
             end do  
            end do
           end do 
         
         end do !m

        do i=1, nmaps
         do j=1, nmaps   
           P%Ps(i,j)%Cl(l,1:numMaps(i),1:numMaps(j)) =Cl(j,i,1:numMaps(i),1:numMaps(j))/(2*l+1)
         end do
        end do 

     end do !l
     deallocate(Cl)

     deallocate(A%alms)
     if (nT/=0) deallocate(AT%alms)
#ifdef MPIPIX
     if (DebugMsgs>0) print *,code // ' Time: ', GeteTime() - iniTime
#endif
    end if
    
  end subroutine maparray2crosspowers   


  subroutine maparray2packedscalalms(H,inlmax, maps, outalm, in_nmap, dofree)
    !=======================================================================
    ! Compute power array 1/(2l+1)\sum_m <alm(1)alm(2)^*> for all maps
    !=======================================================================
    !Allocates outalm%alms
    use MPIStuff
    Type (HealpixInfo) :: H
    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    Type (HealpixMapArray), target :: maps(:)
    Type(HealpixPackedScalAlms) :: outalm
    integer, intent(in) :: in_nmap
    integer nmaps
    logical dofree
    
    COMPLEX(SPC), DIMENSION(:,:), pointer :: alm2
    Type (HealpixMapArray), dimension(:), pointer ::  map2N,map2S
    
    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax, a_ix       ! alm related
    INTEGER(I4B) :: nph, kphi0   ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_0, lam_1, lam_2
    REAL(DP) ::  f2m, corfac

    CHARACTER(LEN=*), PARAMETER :: code = 'MAPARRAY2PACKEDSCALALMS'
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: phas_n
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: phas_s
    INTEGER(I4B) :: mmax_ring, status, nalms
    COMPLEX(DPC), dimension(:), allocatable :: phas_p, phas_m

    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    integer lmin !, par_lm
    integer map_ix
#ifdef MPIPIX
    Type (HealpixMapArray), pointer :: amap
    double precision Initime
#endif
    !=======================================================================


    nsmax = H%nside
    nlmax = inlmax
    nmaps = in_nmap
     
#ifdef MPIPIX
    if (H%MpiId==0) then 
      if(DebugMsgs>1) print *,code //': Sending to farm '
      IniTime = GeteTime()
      call SendMessages(H,code)
    end if

    StartTime = getetime()    
    call SyncInts(nlmax,nmaps)

   if (H%MpiId/=0) then
       allocate(map2N(nmaps),map2S(nmaps))
       do map_ix=1,nmaps
        allocate(map2N(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,1))
        allocate(map2S(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,1))
       end do
   else
     map2N => maps
     map2S=>maps
   end if

   do map_ix=1,nmaps
    if (H%MpiId==0) then
      amap => maps(map_ix)
     else
      amap => map2N(map_ix)
     end if 
    call MPI_SCATTERV(amap%M,H%North_Size, H%North_Start, &
       SP_MPI, map2N(map_ix)%M(H%North_Start(H%MpiId),1),H%North_Size(H%MpiId),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_SCATTERV(amap%M,H%South_Size, H%South_Start, &
       SP_MPI, map2S(map_ix)%M(H%South_Start(H%MpiId),1),H%South_Size(H%MpiId),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
   end do
   if (H%MpiId==0 .and. dofree) then
       allocate(map2N(nmaps),map2S(nmaps))
       do map_ix=1,nmaps
        allocate(map2N(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,1))
        allocate(map2S(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,1))
        map2N(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,1) = &
         maps(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,1)
        map2S(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,1) = &
         maps(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,1) 
        deallocate(maps(map_ix)%M)
       end do
    end if
    if(DebugMsgs>1) print *,code //' Scattered ',H%MpiId, GeteTime() - StartTime
#else
    map2N => maps
    map2S => maps
#endif

    ALLOCATE(phas_n(0:nlmax,nmaps),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_s(0:nlmax,nmaps),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')
  
    ALLOCATE(phas_p(nmaps),phas_m(nmaps)) 
  
    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

    call HealpixInitRecfac(H,nlmax)

    nalms = ((nlmax+1)*(nlmax+2))/2
    
    allocate(outalm%alms(nmaps,nalms), stat=status)
    if (status /= 0) call die_alloc(code,'alm2')
    alm2 => outalm%alms   
  
    alm2 = 0

    omega_pix = pi / (3.0_dp * nsmax * nsmax)

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)

       phas_n = 0._dp 
       phas_s = 0._dp
       
       if (ith .le. nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif

       mmax_ring = get_mmax(nlmax,sth) 

       do map_ix = 1, nmaps
        ring(0:nph-1) = map2N(map_ix)%M(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,1) * H%w8ring_TQU(ith,1)
        call spinring_analysis(H,nlmax, ring, nph, phas_n(:,map_ix), kphi0, mmax_ring)
  
        if (ith .lt. 2*nsmax) then
          ring(0:nph-1) = map2S(map_ix)%M(H%istart_south(ith):H%istart_south(ith)+nph-1,1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_s(:,map_ix), kphi0,mmax_ring)
        endif
       end do

 
       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

          lam_mm = sq4pi_inv * omega_pix ! lambda_00 * norm
          scalem=1 
          a_ix = 0
          do m = 0, mmax_ring
             f2m = 2.0_dp * m
         
             !           ---------- l = m ----------
             if (m .ge. 1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm).lt.UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif
             corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
             
             phas_p=phas_n(m,:) + phas_s(m,:)
             phas_m=phas_n(m,:) - phas_s(m,:)
             
             a_ix = a_ix + 1
             alm2(:,a_ix) = alm2(:,a_ix) + corfac * phas_p
             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec

             lmin = l_min_ylm(m, sth)
             do l = m+1, nlmax-1, 2
   
                a_ix = a_ix+1

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
      
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec
      
                if (l >= lmin) then
                 alm2(:,a_ix) = alm2(:,a_ix) + (lam_1*corfac) * phas_m
                 alm2(:,a_ix+1) = alm2(:,a_ix+1) + (lam_2*corfac) * phas_p
                end if

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_ix = a_ix+1
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_1+lam_2) .gt. OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                elseif (abs(lam_1+lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                endif
             enddo ! loop on l
             if (mod(nlmax-m,2)==1) then
                a_ix = a_ix+1
                alm2(:,a_ix) = alm2(:,a_ix) + (lam_2*corfac) * phas_m
             end if
          enddo ! loop on m
    enddo ! loop on theta


    DEALLOCATE(phas_n, phas_s)
    DEALLOCATE(phas_p, phas_m)
    DEALLOCATE(ring)
    call HealpixFreeRecFac(H)
    
#ifdef MPIPIX
    if (H%MpiId>0 .or. dofree) then
      do map_ix=1,nmaps
       deallocate(map2N(map_ix)%M,map2S(map_ix)%M)
      end do
      if (H%MpiId>0) deallocate(map2S,map2N)
    end if
    
    StartTime = Getetime()
    if (H%MpiId==0) then
      call MPI_REDUCE(MPI_IN_PLACE,alm2,size(alm2),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
      print *,code // ' Time: ', GeteTime() - IniTime
    else
     call MPI_REDUCE(alm2,MPI_IN_PLACE,size(alm2),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
     deallocate(outalm%alms)
    end if
    if (DebugMsgs>1) print *,code //': Time at reduce ', H%MpiId, GeteTime() -StartTime
#else
   if (dofree) then
      do map_ix=1,nmaps
       deallocate(maps(map_ix)%M)
      end do
   end if    
#endif

  END subroutine maparray2packedscalalms


  subroutine maparray2packedscalalms1(H,nlmax, maps, outalm, nmaps)
    !=======================================================================
    ! Compute power array 1/(2l+1)\sum_m <alm(1)alm(2)^*> for all maps
    !=======================================================================
    !Allocates outalm%alms
    !This is non-MPI version
    use MPIStuff
    Type (HealpixInfo) :: H
    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax
    Type (HealpixMapArray), target :: maps(:)
    Type(HealpixPackedScalAlms) :: outalm
    integer, intent(in) :: nmaps
    CHARACTER(LEN=*), PARAMETER :: code = 'MAPARRAY2PACKEDSCALALMS1'
    
    COMPLEX(SPC), DIMENSION(:,:), pointer :: alm2
    
    INTEGER(I4B) :: l, m, ith, scalem, scalel,  a_ix       ! alm related
    INTEGER(I4B) :: nph, kphi0   ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_0, lam_1, lam_2
    REAL(DP) ::  f2m, corfac

    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: phas_n
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: phas_s
    INTEGER(I4B) :: mmax_ring, status, nalms
    COMPLEX(DPC), dimension(:), allocatable :: phas_p, phas_m

    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    integer lmin !, par_lm
    integer map_ix
    integer istart_north, istart_south
    !=======================================================================


    nsmax = H%nside
     
 
    ALLOCATE(phas_n(0:nlmax,nmaps),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_s(0:nlmax,nmaps),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')
  
    ALLOCATE(phas_p(nmaps),phas_m(nmaps)) 
  
    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

    call HealpixInitRecfac(H,nlmax)

    nalms = ((nlmax+1)*(nlmax+2))/2
    
    allocate(outalm%alms(nmaps,nalms), stat=status)
    if (status /= 0) call die_alloc(code,'alm2')
    alm2 => outalm%alms   
  
    alm2 = 0

    istart_north = 0
    istart_south = 12*nsmax**2

    omega_pix = pi / (3.0_dp * nsmax * nsmax)

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    do ith = 1, 2*nsmax

       phas_n = 0._dp 
       phas_s = 0._dp
       
       if (ith .le. nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif

       mmax_ring = get_mmax(nlmax,sth) 

       do map_ix = 1, nmaps
         ring(0:nph-1) = maps(map_ix)%M(istart_north:istart_north+nph-1,1)* H%w8ring_TQU(ith,1)
         call spinring_analysis(H,nlmax, ring, nph, phas_n(:,map_ix), kphi0, mmax_ring)
       end do
       istart_north = istart_north + nph
       istart_south = istart_south - nph

      if (ith .lt. 2*nsmax) then
        do map_ix = 1, nmaps
          ring(0:nph-1) = maps(map_ix)%M(istart_south:istart_south+nph-1,1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_s(:,map_ix), kphi0,mmax_ring)
        end do
       endif
 
 
       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

          lam_mm = sq4pi_inv * omega_pix ! lambda_00 * norm
          scalem=1 
          a_ix = 0
          do m = 0, mmax_ring
             f2m = 2.0_dp * m
         
             !           ---------- l = m ----------
             if (m .ge. 1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm).lt.UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif
             corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
             
             phas_p=phas_n(m,:) + phas_s(m,:)
             phas_m=phas_n(m,:) - phas_s(m,:)
             
             a_ix = a_ix + 1
             alm2(:,a_ix) = alm2(:,a_ix) + corfac * phas_p
             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec

             lmin = l_min_ylm(m, sth)
             do l = m+1, nlmax-1, 2
   
                a_ix = a_ix+1

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
      
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec
      
                if (l >= lmin) then
                 alm2(:,a_ix) = alm2(:,a_ix) + (lam_1*corfac) * phas_m
                 alm2(:,a_ix+1) = alm2(:,a_ix+1) + (lam_2*corfac) * phas_p
                end if

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_ix = a_ix+1
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_1+lam_2) .gt. OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                elseif (abs(lam_1+lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
                endif
             enddo ! loop on l
             if (mod(nlmax-m,2)==1) then
                a_ix = a_ix+1
                alm2(:,a_ix) = alm2(:,a_ix) + (lam_2*corfac) * phas_m
             end if
          enddo ! loop on m
    enddo ! loop on theta


    DEALLOCATE(phas_n, phas_s)
    DEALLOCATE(phas_p, phas_m)
    DEALLOCATE(ring)
    call HealpixFreeRecFac(H)

  END subroutine maparray2packedscalalms1


  subroutine maparray2packedpolalms(h, inlmax, maps,outTEB, in_nmap, dofree)
    use mpistuff
    type (healpixinfo) :: h
    integer(i4b), intent(in) :: inlmax
    Type (HealpixMapArray), target :: maps(:)
    Type(HealpixPackedAlms) :: outTEB
    integer, intent(in) :: in_nmap
    logical, intent(in) :: dofree

    integer(i4b)  :: nsmax, nmaps
    complex(spc), dimension(:,:,:), pointer :: TEB
    Type (HealpixMapArray), dimension(:), pointer ::  map2N,map2S


    integer :: nlmax !base integer type for mpi compat
    integer(i4b) :: l, m, ith, scalem, scalel        ! alm related
    integer(i4b) :: nph, kphi0 ! map related

    real(dp) :: omega_pix
    real(dp) :: cth, sth, dth1, dth2, dst1
    real(dp) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    real(dp) :: fl, fm, f2m, fm2, fm2fac,corfac
    real(dp) :: c_on_s2, fm_on_s2, one_on_s2
    real(dp) :: lambda_w, lambda_x, a_w, b_w, a_x
  
    character(len=*), parameter :: code = 'MAPARRAY2PACKEDPOLALMS'
    complex(dpc), dimension(:,:), allocatable :: phas_nq, phas_nu
    complex(dpc), dimension(:,:), allocatable :: phas_sq, phas_su
    complex(dpc), dimension(:,:), allocatable :: phas_n
    complex(dpc), dimension(:,:), allocatable :: phas_s
    complex(dpc), dimension(:), allocatable :: phas_Qp,phas_Qm,phas_UM,phas_Up,phas_p,phas_m
    complex(dpc), dimension(:), allocatable  :: Iphas_Qp,Iphas_Qm,Iphas_UM,Iphas_Up
    

    integer(i4b) mmax_ring, status, par_lm, a_ix
    integer map_ix, nalms, lmin
    real(dp), dimension(:), allocatable :: lam_fact
    real(dp), dimension(:),   allocatable :: ring
    real(dp), dimension(:),   allocatable :: normal_l
    real(dp), dimension(:), allocatable :: twocthlm1, lfac
    
#ifdef MPIPIX
    Type (HealpixMapArray), pointer :: amap
    integer i
    double precision initime
#endif
    !=======================================================================

    nsmax = h%nside
    nmaps = in_nmap
    nlmax = inlmax
     
#ifdef MPIPIX
     if (h%MpiId==0) then 
      if(debugmsgs>0) print *,code //': sending to farm '
      call sendmessages(h,code)
      map2N => maps
      map2S=>maps
    end if

    starttime = getetime()   
    initime = starttime
    call SyncInts(nlmax,nmaps)
    
    if (H%MpiId /=0) then
       allocate(map2N(nmaps),map2S(nmaps))
       do map_ix=1,nmaps
        allocate(map2N(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,3))
        allocate(map2S(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,3))
       end do
    end if
 
   do map_ix=1,nmaps
     if (H%MpiId==0) then
      amap => maps(map_ix)
     else
      amap => map2N(map_ix)
     end if 
         do i=1,3
          call mpi_scatterv(amap%M(:,i),h%north_size, h%north_start, &
            sp_mpi, map2N(map_ix)%M(h%north_start(h%MpiId),i),h%north_size(h%MpiId),sp_mpi, 0 ,mpi_comm_world, ierr)
          call mpi_scatterv(amap%M(:,i),h%south_size, h%south_start, &
            sp_mpi, map2S(map_ix)%M(h%south_start(h%MpiId),i),h%south_size(h%MpiId),sp_mpi, 0 ,mpi_comm_world, ierr)
         end do
    end do
    if(debugmsgs>1) print *,code //' scattered ',h%MpiId, getetime() - starttime
   if (H%MpiId==0 .and. dofree) then
       allocate(map2N(nmaps),map2S(nmaps))
       do map_ix=1,nmaps
        allocate(map2N(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,3))
        allocate(map2S(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,3))
        do i=1,3
        map2N(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,i) = &
         maps(map_ix)%M(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,i)
        map2S(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,i) = &
         maps(map_ix)%M(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,i) 
        end do   
        deallocate(maps(map_ix)%M)
       end do
   end if
#else
    map2N => maps
    map2S => maps
#endif

    nalms = ((nlmax+1)*(nlmax+2))/2   

    allocate(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    allocate(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')
    allocate(twocthlm1(0:nlmax))
    allocate(lfac(nlmax))
    
    allocate(phas_n(0:nlmax,nmaps), phas_nq(0:nlmax,nmaps), phas_nu(0:nlmax,nmaps)) 
    allocate(phas_s(0:nlmax,nmaps), phas_sq(0:nlmax,nmaps), phas_su(0:nlmax,nmaps)) 
    
    ALLOCATE(phas_Qp(nmaps), phas_Qm(nmaps), phas_UM(nmaps), phas_Up(nmaps), phas_p(nmaps), phas_m(nmaps), &
            Iphas_Qp(nmaps), Iphas_Qm(nmaps), Iphas_UM(nmaps), Iphas_Up(nmaps)) 
  
    allocate(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

   call healpixinitrecfac(h,nlmax)
   call getlamfact(lam_fact, nlmax)
   lam_fact = lam_fact*2

    allocate(outTEB%alms(nmaps,3,nalms), stat=status)
    if (status /= 0) call die_alloc(code,'TEB')
    TEB => outTEB%alms   
    TEB=0
       
    omega_pix = pi / (3 * nsmax * real(nsmax,dp))

    normal_l = 0.0_dp
    do l = 2, nlmax
        fl = dble(l)
        normal_l(l) = eb_sign * sqrt( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo

    dth1 = 1.0_dp / (3.0_dp*dble(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*dble(nsmax))
    dst1 = 1.0_dp / (sqrt(6.0_dp) * dble(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    
    do ith = h%ith_start(h%MpiId), h%ith_end(h%MpiId)

       phas_nq=0; phas_sq=0;phas_nu=0;phas_su=0; phas_n=0; phas_s=0

       if (ith  <=  nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - dble(ith)**2 * dth1
          sth = sin( 2.0_dp * asin( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = mod(ith+1-nsmax,2)
          cth = dble(2*nsmax-ith) * dth2
          sth = dsqrt((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       do l=1,nlmax
        twocthlm1(l) = cth*real(2*(l-1),dp) !! 2(l-1)cos(theta)
        lfac(l) = -real(2*l,dp)*one_on_s2 -real(l,dp)*real(l-1,dp)
       end do

       mmax_ring = get_mmax(nlmax,sth) 

       do map_ix = 1, nmaps 
          ring(0:nph-1) = map2N(map_ix)%M(h%istart_north(ith-1):h%istart_north(ith-1)+nph-1,1) * h%w8ring_TQU(ith,1)
          call spinring_analysis(h,nlmax, ring, nph, phas_n(:,map_ix), kphi0, mmax_ring)
          ring(0:nph-1) = map2N(map_ix)%M(h%istart_north(ith-1):h%istart_north(ith-1)+nph-1,2) * h%w8ring_TQU(ith,2)
          call spinring_analysis(h,nlmax, ring, nph, phas_nq(:,map_ix), kphi0, mmax_ring)
          ring(0:nph-1) = map2N(map_ix)%M(h%istart_north(ith-1):h%istart_north(ith-1)+nph-1,3) * h%w8ring_TQU(ith,3) 
          call spinring_analysis(h,nlmax, ring, nph, phas_nu(:,map_ix), kphi0, mmax_ring)
       
        if (ith  <  2*nsmax ) then
          ring(0:nph-1) = map2S(map_ix)%M(h%istart_south(ith):h%istart_south(ith)+nph-1,1) * h%w8ring_TQU(ith,1)
          call spinring_analysis(h,nlmax, ring, nph, phas_s(:,map_ix), kphi0, mmax_ring)
          ring(0:nph-1) = map2S(map_ix)%M(h%istart_south(ith):h%istart_south(ith)+nph-1,2) * h%w8ring_TQU(ith,2)
          call spinring_analysis(h,nlmax, ring, nph, phas_sq(:,map_ix), kphi0, mmax_ring)
          ring(0:nph-1) = map2S(map_ix)%M(h%istart_south(ith):h%istart_south(ith)+nph-1,3) * h%w8ring_TQU(ith,3)
          call spinring_analysis(h,nlmax, ring, nph, phas_su(:,map_ix), kphi0, mmax_ring)
        endif
       end do
       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

          lam_mm = sq4pi_inv * omega_pix 
          scalem=1
          a_ix = 0
          do m = 0, mmax_ring
             fm  = dble(m)
             f2m = 2.0_dp * fm
             fm2 = fm * fm
             
             fm2fac= 2._dp* fm2 * one_on_s2
             
             fm_on_s2 = fm * one_on_s2

             !           ---------- l = m ----------
             par_lm = 1   ! = (-1)^(l+m+s)
             if (m  >=  1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm) < unflow) then
                lam_mm=lam_mm*ovflow
                scalem=scalem-1
             endif

             a_ix = a_ix+1

             corfac = scalefactor(scalem)*lam_mm/ovflow
             lam_lm = corfac
          
             lmin = max(2,l_min_ylm(m, sth))
             phas_Qp=phas_nQ(m,:) + phas_sQ(m,:)
             phas_Qm=phas_nQ(m,:) - phas_sQ(m,:)
             phas_Up=phas_nU(m,:) + phas_sU(m,:)
             phas_Um=phas_nU(m,:) - phas_sU(m,:)
             Iphas_Qp = (0.0_dp, 1._dp)*phas_Qp
             Iphas_Qm = (0.0_dp, 1._dp)*phas_Qm
             Iphas_Up = (0.0_dp, 1._dp)*phas_Up
             Iphas_Um = (0.0_dp, 1._dp)*phas_Um

             phas_p= phas_n(m,:) + phas_s(m,:)
             phas_m= phas_n(m,:) - phas_s(m,:)


             TEB(:,1,a_ix) = TEB(:, 1,a_ix) + lam_lm * phas_p
             if (m >= lmin) then
            
              lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) *( 2.0_dp * one_on_s2 - 1.0_dp )
              lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
                 TEB(:, 2,a_ix) = TEB(:,2,a_ix) &
                  &                 + lambda_w * phas_Qp + lambda_x*Iphas_Um

                 TEB(:, 3,a_ix) = TEB(:, 3,a_ix) &
                  &                 + lambda_w * phas_Up - lambda_x*Iphas_Qm
            
             end if

             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = h%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec
         
             do l = m+1, nlmax
                par_lm = - par_lm  ! = (-1)^(l+m)
                lam_lm1m=lam_lm ! actual lambda_l-1,m (useful for polarisation)
                lam_lm   = lam_2*corfac ! actual lambda_lm (ovflow factors removed)
                
                a_ix = a_ix + 1

                if (par_lm==1) then
                 TEB(:,1,a_ix) = TEB(:,1,a_ix) + lam_lm * phas_p
                else
                 TEB(:,1,a_ix) = TEB(:,1,a_ix) + lam_lm * phas_m
                end if

             if (l>=lmin .and. corfac /= 0) then
                 !corfac=0 guarantees lam(l-1) is also v close to zero

!                a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)

                 a_w =   fm2fac  + lfac(l)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  twocthlm1(l) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)

       
                if (par_lm==1) then

                 TEB(:,2,a_ix) = TEB(:,2,a_ix) &
                               + lambda_w * phas_qp + lambda_x * Iphas_Um 
                 TEB(:,3,a_ix) = TEB(:,3,a_ix) &
                              +  lambda_w * phas_Up - lambda_x * Iphas_Qm
                
                 else                

                 TEB(:,2,a_ix) = TEB(:,2,a_ix) &
                               + lambda_w * phas_Qm + lambda_x*Iphas_Up
                 TEB(:,3,a_ix) = TEB(:,3,a_ix) &
                              +  lambda_w * phas_Um - lambda_x*Iphas_Qp
       
                 end if
              end if ! l allowed by spin or zero

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_rec = h%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_2)  >  ovflow) then
                   lam_0=lam_0/ovflow
                   lam_1=lam_1/ovflow
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = scalefactor(scalem+scalel)*lam_mm/ovflow
                elseif (abs(lam_2)  <  unflow) then
                   lam_0=lam_0*ovflow
                   lam_1=lam_1*ovflow
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = scalefactor(scalem+scalel)*lam_mm/ovflow
                endif

             enddo ! loop on l
          enddo ! loop on m
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixfreerecfac(h)
    deallocate(lam_fact,lfac)
    deallocate(normal_l,twocthlm1)
    deallocate(phas_nq,phas_nu)
    deallocate(phas_sq,phas_su)
    deallocate(phas_n)
    deallocate(phas_s)
    deallocate(phas_Qp,phas_Qm,phas_UM,phas_Up,phas_p,phas_m,Iphas_Qp,Iphas_Qm,Iphas_UM,Iphas_Up)
    deallocate(ring)
#ifdef MPIPIX
    if (H%MpiId>0 .or. dofree) then
      do map_ix=1,nmaps
       deallocate(map2N(map_ix)%M,map2S(map_ix)%M)
      end do
      if (H%MpiId>0) deallocate(map2S,map2N)
    end if
    starttime = getetime()
    if (H%MpiId==0) then
     call MPI_REDUCE(MPI_IN_PLACE,TEB,size(TEB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    else
     call MPI_REDUCE(TEB,MPI_IN_PLACE,size(TEB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
     deallocate(outTEB%alms)
    end if
    if (debugmsgs>1) print *,code//' time at reduce ', h%MpiId, getetime() -starttime
    if (debugmsgs>0 .and. h%MpiId==0) print *,code //' time: ',getetime() - initime
#else
   if (dofree) then
      do map_ix=1,nmaps
       deallocate(maps(map_ix)%M)
      end do
   end if  

#endif

  end subroutine maparray2packedpolalms



    !=======================================================================
   function scal_at_point(H, cth,sth, phi, alm, nlmax)
    !Get temperature at point from packed alm 
    real scal_at_point   
    Type (HealpixInfo) :: H
    COMPLEX(SPC), INTENT(IN), DIMENSION(:) :: alm
    INTEGER(I4B), INTENT(IN) :: nlmax
    REAL(DP), INTENT(IN) :: cth, sth, phi
    REAL(DP) :: corfac, a_rec, lam_mm, lam_lm, lam_0, f2m, lam_1, lam_2
    INTEGER(I4B) :: par_lm, a_ix, l, m, scalem, scalel          
    integer mmax_ring
    COMPLEX(DPC) :: b_n, factor
    
       lam_mm = sq4pi_inv ! lambda_00
       scalem=1
       a_ix = 0
       mmax_ring = get_mmax(nlmax,sth) 

       do m = 0, mmax_ring
          f2m = 2.0_dp * m

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m)
          if (m >= 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*corfac/OVFLOW
  
          lam_lm = corfac
          a_ix = a_ix + 1
          b_n = lam_lm * alm(a_ix)
          
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp 
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac ! Remove OVFLOW-factors 
             a_ix = a_ix + 1
             factor = lam_lm * alm(a_ix)
             b_n = b_n + factor

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) > OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2) .lt. UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif
          enddo

          if (m==0) then
            scal_at_point = real(b_n)  
          else
            scal_at_point =  scal_at_point + 2*real(b_n*cmplx(cos(m*phi),sin(m*phi)))  
          end if

      end do

    end function scal_at_point


    !=======================================================================
subroutine scalalm2LensedMap(H, inlmax, alm, grad_phi_map, map)
  !Get lensed map by brute force
  !No FFT so does not scale with ln, so slow
  use MPIstuff
  Type (HealpixInfo) :: H
  INTEGER(I4B) :: nsmax
  INTEGER(I4B), INTENT(IN) :: inlmax
  COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
  REAL(SP),     INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map
  COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
  REAL(SP),     DIMENSION(:), pointer :: map2
  COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
  COMPLEX(SPC), DIMENSION(:), allocatable :: alm2
  INTEGER(I4B) :: ith, nlmax          ! alm related
  INTEGER(I4B) :: nph, kphi0                         ! map related
  REAL(DP) :: cth0,sth0,cth, sth, dth1, dth2, dst1
  real(DP) :: phi
  REAL(DP) :: grad_len, sinc_grad_len
  CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2LENSEDMAP'
  INTEGER(I4B) :: mmax_ring 
  integer nalms, ring_ix
#ifdef MPIPIX    
  double precision Initime
  integer status
#endif     
  !=======================================================================
  nsmax = H%nside
  nlmax = inlmax
#ifdef MPIPIX
  StartTime = Getetime()
  iniTime = StartTime
  if (H%MpiId==0) then 
    print *,code //': Sending to farm ' 
    call SendMessages(H,code)
  end if
  call SyncInts(nlmax)
#endif
  nalms = ((nlmax+1)*(nlmax+2))/2   
  allocate(alm2(nalms))
  if (H%MpiId==0) then
    call Alm2PackAlm(alm,alm2,nlmax)
    grad_phi => grad_phi_map
  else
    allocate(grad_phi(0:12*H%nside**2-1))    
  end if
#ifdef MPIPIX 
  call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
  call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
  if(DebugMsgs>1) print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
  allocate(map2(0:12*nsmax**2-1), stat = status)    
  if (status /= 0) call die_alloc(code,'map2')
#else
  map2 => map 
#endif
  map2 = 0
  call HealpixInitRecfac(H,nlmax)
  dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
  dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
  dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )
  write(*,*) "begin lensed map", H%ith_end(H%MpiId) !namikawa
  do ith =H%ith_start(H%MpiId), H%ith_end(H%MpiId)  
    write(*,*) ith !namikawa
    !cos(theta) in the pixelisation scheme
    if (ith.lt.nsmax) then  ! polar cap (north)
      cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
      nph = 4*ith
      kphi0 = 1
      sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
    else                   ! tropical band (north) + equator
      cth0 = DBLE(2*nsmax-ith) * dth2
      nph = 4*nsmax
      kphi0 = MOD(ith+1-nsmax,2)
      sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
    endif
    mmax_ring = get_mmax(nlmax,sth0) 
    do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1
!cth = cos(theta + real(grad_phi(ring_ix)))
!sth = sin(theta + real(grad_phi(ring_ix))) 
!phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph + aimag(grad_phi(ring_ix))/sth0
      phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph
      grad_len = abs(grad_phi(ring_ix))
      if (grad_len>0) then
        sinc_grad_len = sin(grad_len)/grad_len
        cth = cos(grad_len)*cth0 - sinc_grad_len*sth0*real(grad_phi(ring_ix)) 
        sth = sqrt((1._dp-cth)*(1._dp+cth))
        if (sth > 1e-10_dp) then
          phi = phi + asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth ))) 
        end if
      else
        cth=cth0
        sth=sth0
      end if
      if (kphi0==1) phi=phi + pi/nph
      map2(ring_ix) = scal_at_point(H, cth,sth,phi,alm2,nlmax)       
    enddo !ring_ix (phi)
    if (ith < 2*nsmax) then
      do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
        !Think of spherical triangle with points 0, (theta0,phi0), (theta, phi)
        phi = (ring_ix-H%istart_south(ith))*2*pi/nph
        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
          sinc_grad_len = sin(grad_len)/grad_len
          cth = -cos(grad_len)*cth0-sinc_grad_len*sth0*real(grad_phi(ring_ix)) 
          sth = sqrt((1._dp-cth)*(1._dp+cth))
          if (sth > 1e-10_dp) then
            phi = phi +asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth ))) 
          end if
        else
          cth=-cth0
          sth=sth0
        end if
        !cth = cos(pi - theta + real(grad_phi(ring_ix)))
        !sth = sin(pi - theta + real(grad_phi(ring_ix))) 
        if (kphi0==1) phi=phi + pi/nph
        map2(ring_ix) = scal_at_point(H, cth,sth,phi,alm2,nlmax)       
      enddo !ring_ix (phi)
    end if !not middle theta
  enddo    ! loop on cos(theta)
  call healpixFreeRecFac(H)
  deallocate(alm2)
#ifdef MPIPIX
  if(DebugMsgs>1) print *,code //' Gather ',H%MpiId
  StartTime = Getetime()
  call MPI_GATHERV(map2(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),SP_MPI, &
  map,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(map2(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),SP_MPI, &
  map,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
  if (DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
  if (DebugMsgs>0 .and. H%MpiId==0) print *,code //' Time :', GeteTime()-IniTime
  if (H%MpiId/=0) deallocate(grad_phi) 
  deallocate(map2)
#endif
end subroutine scalalm2LensedMap

subroutine alm2Lensedmap(H,inlmax, alm_TEB, grad_phi_map,map_TQU)
  use alm_tools
  use MPIstuff
  Type (HealpixInfo) :: H

  INTEGER(I4B), INTENT(IN) :: inlmax
  integer nsmax
  COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
  REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
  COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB

  REAL(SP),     DIMENSION(:,:), pointer :: map2
  COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map

  COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi

  INTEGER(I4B) ::  l, m, ith, scalem, scalel          ! alm related
  INTEGER(I4B) :: nph, kphi0 ! map related
  REAL(DP) :: cth0,sth0

  REAL(DP) :: cth, sth, dth1, dth2, dst1
  REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
  REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
  REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
  REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
  COMPLEX(DPC) :: zi_lam_x
  COMPLEX(DPC) :: factor_1, factor_2
  COMPLEX(DPC) :: b_n, b_Q, b_U, exp_m_phi, gammfac
  REAL(DP) :: Re, Im, grad_len, sinc_grad_len, phi, gamma 
  INTEGER(I4B) :: ring_ix  

  CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDMAP'
  INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax
  integer, parameter :: spin=2

  REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
  integer a_ix, nalms, lmin
#ifdef MPIPIX
  double precision Initime
  integer i
#endif
  !=======================================================================
  !     --- allocates space for arrays ---

  nsmax = H%nside
  nlmax = inlmax

#ifdef MPIPIX
  StartTime = Getetime()
  iniTime = StartTime
  if (H%MpiId==0) then 
    print *,code //': Sending to farm ' 
    call SendMessages(H,code)
  end if
  call SyncInts(nlmax)
#endif

  nalms = ((nlmax+1)*(nlmax+2))/2   
  allocate(TEB(3,nalms))
  if (H%MpiId==0) then
    call TEB2PackTEB(alm_TEB,TEB,nlmax)
    grad_phi => grad_phi_map
  else
    allocate(grad_phi(0:12*H%nside**2-1))    
  end if

#ifdef MPIPIX
  call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
  call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
  if(DebugMsgs>1) print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
  allocate(map2(0:12*nsmax**2-1,3), stat = status)    
  if (status /= 0) call die_alloc(code,'map2')
#else
  map2 => map_TQU 
#endif

  ALLOCATE(lam_fact(nalms),stat = status)    
  if (status /= 0) call die_alloc(code,'lam_fact')

  ALLOCATE(normal_l(0:nlmax),stat = status)    
  if (status /= 0) call die_alloc(code,'normal_l')

  call HealpixInitRecfac(H,nlmax)
  call GetLamFact(lam_fact, nlmax)
  lam_fact = lam_fact * 2 !HealPix polarization def

  normal_l = 0.0_dp
  do l = 2, nlmax
    fl = DBLE(l)
    normal_l(l) = EB_sign * SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
  end do

  dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
  dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
  dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

  do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1
    if (ith < nsmax) then  ! polar cap (north)
      cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
      nph = 4*ith
      kphi0 = 1
      sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
    else                   ! tropical band (north) + equator
      cth0 = DBLE(2*nsmax-ith) * dth2
      nph = 4*nsmax
      kphi0 = MOD(ith+1-nsmax,2)
      sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
    end if
    mmax_ring = get_mmax(nlmax,sth0) 
    do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1
      phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph
      grad_len = abs(grad_phi(ring_ix))
      if (grad_len>0) then
        sinc_grad_len = sin(grad_len)/grad_len
        cth = cos(grad_len)*cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix))
        sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
        phi = phi  +  asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )))
      else
        cth=cth0
        sth=sth0
      end if
      one_on_s2 = 1.0_dp / sth**2 
      c_on_s2 = cth * one_on_s2
      if (kphi0==1) phi=phi + pi/nph
      lam_mm = sq4pi_inv ! lamda_00
      scalem=1
      a_ix = 0
      do m = 0, mmax_ring
        fm  = DBLE(m)
        f2m = 2.0_dp * fm
        fm2 = fm * fm
        fm_on_s2 = fm * one_on_s2
        ! ---------- l = m ----------
        par_lm = 1  ! = (-1)^(l+m+s)
        if (m  >=  1) then ! lambda_0_0 for m>0
          lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
        end if
        if (abs(lam_mm) < UNFLOW) then
          lam_mm=lam_mm*OVFLOW
          scalem=scalem-1
        end if
        corfac = ScaleFactor(scalem)*lam_mm/OVFLOW 
        ! alm_T * Ylm : Temperature
        lam_lm = corfac    !  actual lambda_mm      
        a_ix = a_ix + 1
        b_n = lam_lm * TEB(1,a_ix)
        !l=m special case
        if (m >=2) then
          lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )
          lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) * 2.0_dp * c_on_s2
          zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
          b_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
          b_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
        else
          b_Q=0
          b_U=0
        end if
        ! ---------- l > m ----------
        lam_0 = 0.0_dp
        lam_1 = 1.0_dp
        scalel=0
        a_rec = H%recfac(a_ix)
        lam_2 = cth * lam_1 * a_rec
        lmin = max(2,l_min_ylm(m, sth))
        do l = m+1, nlmax
          par_lm = - par_lm  ! = (-1)^(l+m+s)
          lam_lm1m=lam_lm  ! actual lambda_l-1,m 
          lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
          fl  = DBLE(l)
          fl2 = fl * fl
          a_ix = a_ix + 1
          b_n = b_n + lam_lm * TEB(1,a_ix)
          if (l>=lmin .and. corfac /= 0) then
            a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
            b_w =  c_on_s2 * lam_fact(a_ix)
            a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
            lambda_w = normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
            lambda_x = normal_l(l) * fm_on_s2*( lam_fact(a_ix) * lam_lm1m - a_x)
            zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
            factor_1 =  lambda_w * TEB(2,a_ix)
            factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
            b_Q = b_Q +           factor_1 + factor_2
            factor_1 =   lambda_w * TEB(3,a_ix) 
            factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
            b_U = b_U +           factor_1 - factor_2
          end if
          lam_0 = lam_1 / a_rec
          lam_1 = lam_2
          a_rec = H%recfac(a_ix)
          lam_2 = (cth * lam_1 - lam_0) * a_rec
          if (abs(lam_2)  >  OVFLOW) then
            lam_0=lam_0/OVFLOW
            lam_1=lam_1/OVFLOW
            lam_2 = (cth * lam_1 - lam_0) * a_rec
            scalel=scalel+1
            corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
          else if (abs(lam_2)  <  UNFLOW) then
            lam_0=lam_0*OVFLOW
            lam_1=lam_1*OVFLOW
            lam_2 = (cth * lam_1 - lam_0) * a_rec 
            scalel=scalel-1
            corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
          end if
        end do
        if (m==0) then
          map2(ring_ix,1) =  real(b_n)
          map2(ring_ix,2) =  real(b_Q)
          map2(ring_ix,3) =  real(b_U)  
        else
          exp_m_phi = cmplx(cos(m*phi),sin(m*phi)) 
          map2(ring_ix,1) =  map2(ring_ix,1) + 2*real(b_n*exp_m_phi) 
          map2(ring_ix,2) =  map2(ring_ix,2) + 2*real(b_Q*exp_m_phi)
          map2(ring_ix,3) =  map2(ring_ix,3) + 2*real(b_U*exp_m_phi)  
        end if
      end do
      !Put in factor for change of co-ordinate axes between deflected and original point
      if (grad_len > 1e-20_dp) then
        Re = real(grad_phi(ring_ix))
        Im = aimag(grad_phi(ring_ix))
        gamma = grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
        if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
        gamma = Im/gamma
        !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
        gammfac = cmplx(map2(ring_ix,2),map2(ring_ix,3))* &
        cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
          2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))
        map2(ring_ix,2) = real(gammfac)
        map2(ring_ix,3) = aimag(gammfac)
      end if
    end do !ring_ix (phi)
    if (ith < 2*nsmax) then
      do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
        phi = (ring_ix-H%istart_south(ith))*2*pi/nph
        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
          sinc_grad_len = sin(grad_len)/grad_len
          cth = -cos(grad_len)*cth0 - sinc_grad_len*sth0*real(grad_phi(ring_ix))
          sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
          phi = phi + asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )))
        else
          cth=-cth0
          sth=sth0
        end if
        if (kphi0==1) phi=phi + pi/nph
        one_on_s2 = 1.0_dp / sth**2 
        c_on_s2 = cth * one_on_s2
        lam_mm = sq4pi_inv ! lamda_00
        scalem=1
        a_ix = 0
        do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2
          ! ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
            lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          end if
          if (abs(lam_mm) < UNFLOW) then
            lam_mm=lam_mm*OVFLOW
            scalem=scalem-1
          end if
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
          ! alm_T * Ylm : Temperature
          lam_lm = corfac     !  actual lambda_mm      
          a_ix = a_ix + 1
          b_n = lam_lm * TEB(1,a_ix)
          !l=m special case
          if (m >=2) then
            lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
            lambda_x = ( normal_l(m)*lam_lm * (fm - fm2) ) * 2.0_dp * c_on_s2
            zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
            b_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
            b_U =  lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
          else
            b_Q=0
            b_U=0
          end if
          ! ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          lmin = max(2,l_min_ylm(m, sth))
          do l = m+1, nlmax
            par_lm = - par_lm  ! = (-1)^(l+m+s)
            lam_lm1m=lam_lm  ! actual lambda_l-1,m 
            lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
            fl  = DBLE(l)
            fl2 = fl * fl
            a_ix = a_ix + 1
            b_n = b_n + lam_lm * TEB(1,a_ix)
            if (l>=lmin .and. corfac /= 0) then
              a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
              b_w =  c_on_s2 * lam_fact(a_ix)
              a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
              lambda_w = normal_l(l)*( a_w * lam_lm + b_w * lam_lm1m ) 
              lambda_x = normal_l(l)*fm_on_s2*( lam_fact(a_ix) * lam_lm1m - a_x)
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
              factor_1 =  lambda_w * TEB(2,a_ix)
              factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
              b_Q = b_Q +           factor_1 + factor_2
              factor_1 =   lambda_w * TEB(3,a_ix) 
              factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
              b_U = b_U +           factor_1 - factor_2
            end if
            lam_0 = lam_1 / a_rec
            lam_1 = lam_2
            a_rec = H%recfac(a_ix)
            lam_2 = (cth * lam_1 - lam_0) * a_rec
            if (abs(lam_2)  >  OVFLOW) then
              lam_0=lam_0/OVFLOW
              lam_1=lam_1/OVFLOW
              lam_2 = (cth * lam_1 - lam_0) * a_rec
              scalel=scalel+1
              corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
            else if (abs(lam_2)  <  UNFLOW) then
              lam_0=lam_0*OVFLOW
              lam_1=lam_1*OVFLOW
              lam_2 = (cth * lam_1 - lam_0) * a_rec 
              scalel=scalel-1
              corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
            end if
          end do
          if (m==0) then
            map2(ring_ix,1) =  real(b_n)
            map2(ring_ix,2) =  real(b_Q)
            map2(ring_ix,3) =  real(b_U)  
          else
            exp_m_phi = cmplx(cos(m*phi),sin(m*phi)) 
            map2(ring_ix,1) =  map2(ring_ix,1) + 2*real(b_n*exp_m_phi) 
            map2(ring_ix,2) =  map2(ring_ix,2) + 2*real(b_Q*exp_m_phi)
            map2(ring_ix,3) =  map2(ring_ix,3) + 2*real(b_U*exp_m_phi)  
          end if
        end do
        !Put in factor for change of co-ordinate axes between deflected and original point
        if (grad_len > 1e-20_dp) then
          Re = real(grad_phi(ring_ix))
          Im = aimag(grad_phi(ring_ix))
          gamma = -grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
          !AL: Oct 07, sign of cth0
          if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
          gamma = Im/gamma
          !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
          gammfac = cmplx(map2(ring_ix,2),map2(ring_ix,3))* &
          cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
            2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))
          map2(ring_ix,2) = real(gammfac)
          map2(ring_ix,3) = aimag(gammfac)
        end if
      end do
    end if  
  end do    ! loop on cos(theta)

  !     --------------------
  !     free memory and exit
  !     --------------------
  call HealpixFreeRecfac(H)
  DEALLOCATE(lam_fact)
  DEALLOCATE(normal_l)
  deallocate(TEB)
#ifdef MPIPIX
  if(DebugMsgs>1) print *,code//' Gather ',H%MpiId
  StartTime = Getetime()
  do i=1,3
   call MPI_GATHERV(map2(H%North_Start(H%MpiId),i),H%North_Size(H%MpiId),SP_MPI, &
     map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
   call MPI_GATHERV(map2(H%South_Start(H%MpiId),i),H%South_Size(H%MpiId),SP_MPI, &
     map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
  end do
  if (H%MpiId/=0) deallocate(grad_phi)
  if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
  if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
  deallocate(map2) 
#endif

end subroutine alm2Lensedmap


 !=======================================================================
  subroutine ang2pix_ring8(nside, costheta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map 
    !     resolution parameter nside
    !=======================================================================
    !AL: Uses full I8B integer pixel range, takes in cos(theta) rather than theta
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I8B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  costheta, phi

    INTEGER(KIND=I8B) ::  nl4, jp, jm
    REAL(KIND=DP) ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=I8B) ::  ir, ip, kshift

    z = costheta
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)

    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_dp+tt)
       temp2 = nside*.75_dp*z
       jp = int(temp1-temp2) ! index of  ascending edge line 
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2_I8B) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*int(nside-1,I8B) + nl4*(ir-1) + ip 

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*int(nside,I8B)**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring8

subroutine scalalm2LensedmapInterp(H,inlmax,alm,grad_phi_map,map,nside_factor)
  !AL: Added Oct 2007 
  !Temperature-only lensing by naive pixel remapping to higher-res map, no interpolation 
  use MPIstuff
  Type (HealpixInfo) :: H, H_res
  INTEGER(I4B), INTENT(IN) :: inlmax
  INTEGER(I4B), INTENT(IN), optional :: nside_factor
  integer nsmax
  integer  :: nside_fac = 8
  COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
  COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
  REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map
  REAL(SP), DIMENSION(:), pointer :: map2N, map2S
  REAL(SP), DIMENSION(:), pointer :: high_resN,high_resS
  COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
  COMPLEX(SPC), DIMENSION(:), allocatable :: alm2
  INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
  INTEGER(I4B) :: nph, kphi0 ! map related
  REAL(DP) :: cth, sth, dth1, dth2, dst1
  REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
  REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
  COMPLEX(DPC) :: factor
  COMPLEX(DPC) :: b_n, b_s
  CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2LENSEDMAPINTERP'
  COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
  INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax
  REAL(SP), DIMENSION(:), allocatable ::  ring
  integer(I4B) a_ix, nalms, border_inc, high_th_start,high_th_end, high_nside
  integer(I8B) ipix, ring_ix, resN_start, resS_start
  REAL(DP) :: grad_len, sinc_grad_len, phi 
  REAL(DP) :: cth0, sth0 
#ifdef MPIPIX    
  double precision Initime
#endif !=======================================================================

  !   --- allocates space for arrays ---
  nsmax = H%nside
  nlmax = inlmax
  if (present(nside_factor)) nside_fac = nside_factor
#ifdef MPIPIX
  StartTime = Getetime()
  iniTime = StartTime
  if (H%MpiId==0) then 
    print *,code //': Sending to farm ' 
    call SendMessages(H,code)
  end if
#endif
  if (H%MpiId ==0) then
    high_nside = nsmax*nside_fac
!border_inc=number of high-res pixels we need to go outside zero-lensing border
!1.18 is approx 3Pi/8 which is the large-n_side vertical overdensity of pixels near the 
!equator relative to the average. Could speed by putting in position dependent
!factor. Corrected AL: 30 Sept 2004
    border_inc = int(maxval(abs(real(grad_phi_map)))/PI*4*high_nside*1.18) + 1
  end if
  call SyncInts(nlmax,nside_fac,border_inc)
  nalms = ((nlmax+1)*(nlmax+2))/2   
  allocate(alm2(nalms))
  if (H%MpiId==0) call Alm2PackAlm(alm,alm2,nlmax)
#ifdef MPIPIX 
  call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
  if(DebugMsgs>1) then
    print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
    StartTime = geteTime()
  end if
#endif
  high_nside = H%nside*nside_fac
  call HealpixInitTrig(H_res,high_nside,nlmax)
  H_res%MpiId = 1
  high_th_start = max((H%ith_start(H%MpiId)-1) * nside_fac - border_inc, 1)
  high_th_end = min(H%ith_end(H%MpiId)  * nside_fac + border_inc, 2*high_nside)
  if (high_th_end  < high_nside) then  ! polar cap (north)
    nph = 4*high_th_end
  else                   
    nph = 4*high_nside
  endif
  resN_start = H_res%istart_north(high_th_start-1) 
  allocate(high_resN(0:H_res%istart_north(high_th_end-1)+nph-1 - resN_start))
  if (high_th_start  < high_nside) then  ! polar cap (north)
    nph = 4*high_th_start
  else                   
    nph = 4*high_nside
  endif
  resS_start = H_res%istart_south(min(high_th_end,2*high_nside-1)) 
  allocate(high_resS(0:H_res%istart_south(high_th_start)+nph-1 - resS_start))

  ALLOCATE(ring(0:4*H%nside*nside_fac-1), stat = status)
  if (status /= 0) call die_alloc(code,'ring')
  ! ------------ initiate arrays ----------------
  call HealpixInitRecfac(H,nlmax)
  dth1 = 1.0_dp / (3.0_dp*DBLE(high_nside)**2)
  dth2 = 2.0_dp / (3.0_dp*DBLE(high_nside))
  dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(high_nside) )
  do ith = high_th_start, high_th_end      ! 0 <= cos theta < 1
    !cos(theta) in the pixelisation scheme
    if (ith < high_nside) then  ! polar cap (north)
      cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
      nph = 4*ith
      kphi0 = 1
      sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
    else                   ! tropical band (north) + equator
      cth = DBLE(2*high_nside-ith) * dth2 !cos theta
      nph = 4*high_nside
      kphi0 = MOD(ith+1-high_nside,2)
      sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
    endif
    !        -----------------------------------------------------
    !        for each theta, and each m, computes
    !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
    !        ------------------------------------------------------
    !        lambda_mm tends to go down when m increases (risk of underflow)
    !        lambda_lm tends to go up   when l increases (risk of overflow)
    lam_mm = sq4pi_inv ! lamda_00
    scalem=1
    mmax_ring = get_mmax(nlmax,sth) 
    a_ix = 0
    do m = 0, mmax_ring
      fm  = DBLE(m)
      f2m = 2.0_dp * fm
      fm2 = fm * fm
      !           ---------- l = m ----------
      par_lm = 1  ! = (-1)^(l+m+s)
      if (m  >=  1) then ! lambda_0_0 for m>0
        lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
      endif
      if (abs(lam_mm) < UNFLOW) then
        lam_mm=lam_mm*OVFLOW
        scalem=scalem-1
      endif
      corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
      ! alm_T * Ylm : Temperature
      lam_lm = corfac    !  actual lambda_mm      
      a_ix = a_ix + 1
      b_n = lam_lm * alm2(a_ix)
      b_s = b_n
      !           ---------- l > m ----------
      lam_0 = 0.0_dp
      lam_1 = 1.0_dp
      scalel=0
      a_rec = H%recfac(a_ix)
      lam_2 = cth * lam_1 * a_rec
      do l = m+1, nlmax
        par_lm = - par_lm  ! = (-1)^(l+m+s)
        lam_lm1m=lam_lm  ! actual lambda_l-1,m 
        lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
        fl  = DBLE(l)
        fl2 = fl * fl
        a_ix = a_ix + 1
        factor = lam_lm * alm2(a_ix)
        b_n = b_n + factor
        b_s = b_s + par_lm * factor
        lam_0 = lam_1 / a_rec
        lam_1 = lam_2
        a_rec = H%recfac(a_ix)
        lam_2 = (cth * lam_1 - lam_0) * a_rec
        if (abs(lam_2)  >  OVFLOW) then
          lam_0=lam_0/OVFLOW
          lam_1=lam_1/OVFLOW
          lam_2 = (cth * lam_1 - lam_0) * a_rec
          scalel=scalel+1
          corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
        else if (abs(lam_2)  <  UNFLOW) then
          lam_0=lam_0*OVFLOW
          lam_1=lam_1*OVFLOW
          lam_2 = (cth * lam_1 - lam_0) * a_rec 
          scalel=scalel-1
          corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
        end if
      end do
      b_north(m) = b_n
      b_south(m) = b_s
    end do
    call spinring_synthesis(H_res,nlmax,b_north,nph,ring,kphi0,mmax_ring)   
    high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start) = ring(0:nph-1)
    if (ith  <  2*high_nside) then
      call spinring_synthesis(H_res,nlmax, b_south, nph, ring, kphi0,mmax_ring)
      high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start) = ring(0:nph-1)
    end if
  end do    ! loop on cos(theta)
  !     --------------------
  !     free memory and exit
  !     --------------------
  DEALLOCATE(ring)
  deallocate(alm2)
  deallocate(H_res%trig)
  deallocate(H%recfac,stat= status)
  nullify(H%recfac)
  if (H%MpiId==0) then
    grad_phi => grad_phi_map
  else
    allocate(grad_phi(0:12*H%nside**2-1))    
    grad_phi = 0
  end if
#ifdef MPIPIX
  call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
  if(DebugMsgs>1) then
    print *,code //': Got grad_phi ',H%MpiId, GeteTime() - StartTime
  end if
#endif
#ifdef MPIPIX
  allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1), stat = status)
  if (status /= 0) call die_alloc(code,'map2N')
  allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1),stat = status)
  if (status /= 0) call die_alloc(code,'map2S')
#else
  map2N => map
  map2S => map
#endif
  dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
  dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
  dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )
  do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1
    if (ith < nsmax) then  ! polar cap (north)
      cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
      nph = 4*ith
      kphi0 = 1
      sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
    else                   ! tropical band (north) + equator
      cth0 = DBLE(2*nsmax-ith) * dth2
      nph = 4*nsmax
      kphi0 = MOD(ith+1-nsmax,2)
      sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
    endif
    do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1
      phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph
      grad_len = abs(grad_phi(ring_ix))
      if (grad_len>0) then
        sinc_grad_len = sin(grad_len)/grad_len
        cth = cos(grad_len)*cth0 - sinc_grad_len*sth0*real(grad_phi(ring_ix))
        sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
        phi = phi  + asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )))
      else
        cth=cth0
        sth=sth0
      end if
      if (kphi0==1) phi=phi + pi/nph
      call ang2pix_ring8(high_nside, cth, phi, ipix)
      if (ipix >= H_res%istart_south(2*high_nside-1)) then
        map2N(ring_ix) = high_resS(ipix-resS_start)  
      else
        map2N(ring_ix) = high_resN(ipix-resN_start)
      end if         
    end do !ring_ix
    if (ith < 2*nsmax) then
      do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
        phi = (ring_ix-H%istart_south(ith))*2*pi/nph
        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
          sinc_grad_len = sin(grad_len)/grad_len
          cth =  -cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
          sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
          phi = phi + asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )))
        else
          cth=-cth0
          sth=sth0
        end if
        if (kphi0==1) phi=phi + pi/nph
        call ang2pix_ring8(high_nside, cth, phi, ipix)
        if (ipix >= H_res%istart_south(2*high_nside-1)) then
          map2S(ring_ix) = high_resS(ipix-resS_start)  
        else
          map2S(ring_ix) = high_resN(ipix-resN_Start)
        end if         
      end do
    end if      
  end do !ith
  deallocate(high_resS,high_resN)
  call HealpixFree(H_res)
#ifdef MPIPIX
  if (H%MpiId/=0) then
    deallocate(grad_phi)    
  end if
  if(DebugMsgs>1) print *,code//' Gather ',H%MpiId, GeteTime()-StartTime
  StartTime = Getetime()
  call MPI_GATHERV(map2N(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),SP_MPI, &
  map(:),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
  call MPI_GATHERV(map2S(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),SP_MPI, &
       map(:),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
  if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
  if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
  deallocate(map2N,map2S) 
#endif
end subroutine scalalm2LensedmapInterp


  subroutine alm2LensedmapInterp(H,inlmax, alm_TEB, grad_phi_map, map_TQU, nside_factor)
  !This routine is designed to be used over a cluster of say 30+ nodes, at least
  !for high resolution maps. Without a cluster memory use will probably be > 2GB.
  !It currently just re-maps pixels, without any interpolation or accounting for pixel shape etc
  !nside_fac=4 at nside=1024 is probably sufficient before Planck, nside_fac=8 for Planck at 0.5%.
    use MPIstuff
    Type (HealpixInfo) :: H, H_res

    INTEGER(I4B), INTENT(IN) :: inlmax
    INTEGER(I4B), INTENT(IN), optional :: nside_factor

    integer nsmax
    integer  :: nside_fac = 8
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    REAL(SP), DIMENSION(:,:), pointer :: map2N, map2S
    REAL(SP), DIMENSION(:,:), pointer :: high_resN,high_resS
    COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n, b_s

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDMAPINTERP'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:), allocatable ::  ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer(I4B) a_ix, nalms, border_inc, high_th_start,high_th_end, high_nside
    integer(I8B) ipix, ring_ix, resN_start, resS_start
    COMPLEX(DPC) :: gammfac
    REAL(DP) :: grad_len, sinc_grad_len, phi, gamma  
    REAL(DP) :: Re, Im, cth0, sth0 
#ifdef MPIPIX
    double precision Initime
    integer i
#endif
    !=======================================================================

    !     --- allocates space for arrays ---

     nsmax = H%nside
     nlmax = inlmax
     if (present(nside_factor)) nside_fac = nside_factor

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
#endif

    if (H%MpiId ==0) then
     high_nside = nsmax*nside_fac
 !border_inc=number of high-res pixels we need to go outside zero-lensing border
 !1.18 is approx 3Pi/8 which is the large-n_side vertical overdensity of pixels near the 
 !equator relative to the average. Could speed by putting in position dependent
 !factor. Corrected AL: 30 Sept 2004
     border_inc = int(maxval(abs(real(grad_phi_map)))/ PI *4*high_nside * 1.18) + 1  
    end if

     call SyncInts(nlmax,nside_fac,border_inc)
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(TEB(3,nalms))
     if (H%MpiId==0) then
        call TEB2PackTEB(alm_TEB,TEB,nlmax)
     end if
     
#ifdef MPIPIX
     call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) then
      print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
      StartTime = geteTime()
     end if
#endif

    high_nside = H%nside*nside_fac
    call HealpixInitTrig(H_res,high_nside,nlmax)
    H_res%MpiId = 1
 
    high_th_start = max((H%ith_start(H%MpiId)-1) * nside_fac - border_inc, 1)
    high_th_end  =  min(H%ith_end(H%MpiId)  * nside_fac + border_inc, 2*high_nside)
    if (high_th_end  < high_nside) then  ! polar cap (north)
              nph = 4*high_th_end
           else                   
              nph = 4*high_nside
    endif
    resN_start = H_res%istart_north(high_th_start-1) 
    allocate(high_resN(0:H_res%istart_north(high_th_end-1)+nph-1 - resN_start,3))
    if (high_th_start  < high_nside) then  ! polar cap (north)
              nph = 4*high_th_start
           else                   
              nph = 4*high_nside
    endif
    resS_start = H_res%istart_south(min(high_th_end,2*high_nside-1)) 
    allocate(high_resS(0:H_res%istart_south(high_th_start)+nph-1 - resS_start,3))    

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')

    ALLOCATE(ring(0:4*H%nside*nside_fac-1), stat = status)
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   lam_fact = lam_fact * 2 !HealPix polarization def

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
       normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(high_nside)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(high_nside))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(high_nside) )

    do ith = high_th_start, high_th_end      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
    
       if (ith < high_nside) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*high_nside-ith) * dth2 !cos theta
          nph = 4*high_nside
          kphi0 = MOD(ith+1-high_nside,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       mmax_ring = get_mmax(nlmax,sth) 

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac     !  actual lambda_mm      

          a_ix = a_ix + 1

          b_n = lam_lm * TEB(1,a_ix)
          b_s = b_n

          !l=m special case
          if (m >=2) then
              lambda_w = - 2.0_dp *(normal_l(m) * lam_lm * (fm - fm2) ) * ( one_on_s2 - 0.5_dp )
              lambda_x =  ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              b_n_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_s_Q =  par_lm*(lambda_w * TEB(2,a_ix) - zi_lam_x * TEB(3,a_ix))

              b_n_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
              b_s_U = par_lm*(lambda_w * TEB(3,a_ix) + zi_lam_x * TEB(2,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1

             factor = lam_lm * TEB(1,a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U
          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H_res,nlmax,b_north,nph,ring,kphi0,mmax_ring)   
       high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start,1) = ring(0:nph-1)
       call spinring_synthesis(H_res,nlmax, b_north_Q, nph, ring, kphi0,mmax_ring)
       high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start,2) = ring(0:nph-1)
       call spinring_synthesis(H_res,nlmax, b_north_U, nph, ring, kphi0,mmax_ring)
       high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start,3) = ring(0:nph-1)
  
       if (ith  <  2*high_nside) then
          call spinring_synthesis(H_res,nlmax, b_south, nph, ring, kphi0,mmax_ring)
          high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start,1) = ring(0:nph-1)
          call spinring_synthesis(H_res,nlmax, b_south_Q, nph, ring, kphi0,mmax_ring)
          high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start,2) = ring(0:nph-1)
          call spinring_synthesis(H_res,nlmax, b_south_U, nph, ring, kphi0,mmax_ring)
          high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start,3) = ring(0:nph-1)
       endif

    enddo    ! loop on cos(theta)


    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l, ring)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)

    deallocate(TEB)
    deallocate(H_res%trig)

    deallocate(H%recfac,stat= status)
    nullify(H%recfac)

     if (H%MpiId==0) then
       grad_phi => grad_phi_map
     else
        allocate(grad_phi(0:12*H%nside**2-1))    
        grad_phi = 0
     end if
#ifdef MPIPIX
     call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) then
      print *,code //': Got grad_phi ',H%MpiId, GeteTime() - StartTime
     end if
#endif


#ifdef MPIPIX
     allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,3), stat = status)
     if (status /= 0) call die_alloc(code,'map2N')
     allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,3),stat = status)
     if (status /= 0) call die_alloc(code,'map2S')
#else
     map2N => map_TQU
     map2S => map_TQU
#endif


    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1

       if (ith < nsmax) then  ! polar cap (north)
          cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth0 = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
       endif

       do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1

        phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph

        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
            sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
            phi = phi  + asin(max(-1._dp,min(1._dp,aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )))
        else
         cth=cth0
         sth=sth0
        endif

        if (kphi0==1) phi=phi + pi/nph

        call ang2pix_ring8(high_nside, cth, phi, ipix)
        if (ipix >= H_res%istart_south(2*high_nside-1)) then
          map2N(ring_ix,:) = high_resS(ipix-resS_start,:)  
        else
          map2N(ring_ix,:) = high_resN(ipix-resN_start,:)
        end if         

     !Put in factor for change of co-ordinate axes between deflected and original point
          if (grad_len > 1e-20_dp) then
             Re = real(grad_phi(ring_ix))
             Im = aimag(grad_phi(ring_ix))
             gamma = grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             gammfac = cmplx(map2N(ring_ix,2),map2N(ring_ix,3))* &
               cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
             map2N(ring_ix,2) = real(gammfac)
             map2N(ring_ix,3) = aimag(gammfac)
         end if

       end do !ring_ix
      
      if (ith < 2*nsmax) then
       do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
            phi = (ring_ix-H%istart_south(ith))*2*pi/nph
            grad_len = abs(grad_phi(ring_ix))
            if (grad_len>0) then
                sinc_grad_len = sin(grad_len)/grad_len
                cth =  -cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
                sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
                phi = phi +  asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )
            else
             cth=-cth0
             sth=sth0
            endif
           if (kphi0==1) phi=phi + pi/nph
           call ang2pix_ring8(high_nside, cth, phi, ipix)
           if (ipix >= H_res%istart_south(2*high_nside-1)) then
              map2S(ring_ix,:) = high_resS(ipix-resS_start,:)  
           else
              map2S(ring_ix,:) = high_resN(ipix-resN_Start,:)
           end if         

          
          if (grad_len > 1e-20_dp) then
             Re = real(grad_phi(ring_ix))
             Im = aimag(grad_phi(ring_ix))
             gamma = -grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
              !AL: Oct 07: fixed sign of cth0 for southern
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             gammfac = cmplx(map2S(ring_ix,2),map2S(ring_ix,3))* &
               cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
             map2S(ring_ix,2) = real(gammfac)
             map2S(ring_ix,3) = aimag(gammfac)
         end if

       end do
      end if      

    end do !ith

    deallocate(high_resS,high_resN)
    call HealpixFree(H_res)

#ifdef MPIPIX
    if (H%MpiId/=0) then
        deallocate(grad_phi)    
    end if
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiId, GeteTime()-StartTime
    StartTime = Getetime()
    do i=1,3
    call MPI_GATHERV(map2N(H%North_Start(H%MpiId),i),H%North_Size(H%MpiId),SP_MPI, &
       map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2S(H%South_Start(H%MpiId),i),H%South_Size(H%MpiId),SP_MPI, &
       map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end do
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2N,map2S) 
#endif

  end subroutine alm2LensedmapInterp



    subroutine cyl_interp_init(H, Grad, grad_phi_map,nside_fac, high_nside, n_phi_cyl, dphi_cyl,dtheta_cyl,&
         cyl_start_ix, cyl_end_ix)
     use MPIstuff
     Type (HealpixInfo) :: H  
     Type (LensGradients) :: Grad
     COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
     real, intent(in) :: nside_fac
     integer, intent(out) :: n_phi_cyl, high_nside
     integer, intent(out) :: cyl_start_ix, cyl_end_ix
     real(dp), intent(out) :: dphi_cyl,  dtheta_cyl
     real(dp) :: dth1, dth2, cth, cth2, sth
     integer(I4B) :: border_N, border_S  
     integer(I4B), parameter :: phi_extra_fac = 1
     integer ith, nph
     CHARACTER(LEN=*), PARAMETER :: code = 'cyl_interp_init'
#ifdef MPIPIX
     integer status
#endif

     high_nside = max(H%nside,8*nint((H%nside*nside_fac)/8))

#ifdef MPIPIX
    if(DebugMsgs>1 .and. H%MpiId==0) print *, code //' interpolation with nside = ', high_nside 
    allocate(Grad%grad_phiN(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1), stat = status)
    if (status /= 0) call die_alloc(code,'grad_phiN')
    allocate(Grad%grad_phiS(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1),stat = status)
     if (status /= 0) call die_alloc(code,'grad_phiS')

    call MPI_SCATTERV(grad_phi_map,H%North_Size, H%North_Start, &
       CSP_MPI, Grad%grad_phiN(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)

    call MPI_SCATTERV(grad_phi_map,H%South_Size, H%South_Start, &
       CSP_MPI, Grad%grad_phiS(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *, code //' Scattered grad_phi',H%MpiId, GeteTime() - StartTime

#else
     Grad%grad_phiN => grad_phi_map
     Grad%grad_phiS => grad_phi_map
#endif
 
 
     
!Get angular range of this chunk to split up
     dth1 = 1.0_dp / (3.0_dp*DBLE(H%nside)**2)
     dth2 = 2.0_dp / (3.0_dp*DBLE(H%nside))
    
     if (H%ith_start(H%MpiId) < H%nside) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(H%ith_start(H%MpiId) )**2 * dth1  !cos theta
     else                   ! tropical band (north) + equator
          cth = DBLE(2*H%nside-H%ith_start(H%MpiId) ) * dth2 !cos theta
     endif
    
     if (H%ith_end(H%MpiId) < H%nside) then  ! polar cap (north)
          cth2 = 1.0_dp  - DBLE(H%ith_end(H%MpiId) )**2 * dth1  !cos theta
     else                   ! tropical band (north) + equator
          cth2 = DBLE(2*H%nside-H%ith_end(H%MpiId) ) * dth2 !cos theta                              
     endif
     sth = sqrt((1._dp-cth2)*(1._dp+cth2))
 
     dtheta_cyl = PI/real(2._dp*high_nside,DP)
 
     n_phi_cyl = max(128,NearestFastFFTnum(nint(4*high_nside*sth*phi_extra_fac)))
    
    dphi_cyl = 2._dp*pi/real(n_phi_cyl,DP)
!first pixel centre is at dtheta_cyl/2 
!equicylindric pixels numbered in theta from 0 -> high_nside -1  
 
!Find the borders in pixels need to go top and bottom of strip because of deflections out of strip
    border_N=0
    border_S=0
    do ith = H%ith_start(H%MpiId), min(2*H%nside-1,H%ith_end(H%MpiId))      !
       if (ith < H%nside) then  ! polar cap (north)
          nph = 4*ith
       else                   ! tropical band (north) + equator
          nph = 4*H%nside
       endif
         
       border_N = max(border_N,  &
          1-int(minval(real(Grad%grad_phiN(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1))) &
                                / PI *2*high_nside  + (ith - H%ith_start(H%MpiId))*nside_fac/2.5)  )

       border_N = max(border_N,&
          1+int(maxval(real(Grad%grad_phiS(H%istart_South(ith):H%istart_south(ith)+nph-1))) &
                                / PI *2*high_nside  - (ith - H%ith_start(H%MpiId))*nside_fac/2.5) )

    
       border_S = max(border_S, &
          1+int(maxval(real(Grad%grad_phiN(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1))) &
                                / PI *2*high_nside  - (H%ith_end(H%MpiId)-ith)*nside_fac/2.5) )
    
       border_S = max(border_S, &
          1-int(minval(real(Grad%grad_phiS(H%istart_South(ith):H%istart_south(ith)+nph-1))) &
                                / PI *2*high_nside  + (H%ith_end(H%MpiId)-ith)*nside_fac/2.5) )
    end do                             

#ifdef MPIPIX
   if (DebugMsgs > 1) print *, H%MpiId, 'borders = ', border_S, border_N
#endif
    cyl_start_ix =  int(dacos(cth) / dtheta_cyl) - border_N
    cyl_end_ix   =  int(dacos(cth2) / dtheta_cyl+1)+  border_S
        
    end subroutine cyl_interp_init


   subroutine scalalm2LensedmapInterpCyl(H,inlmax, alm, grad_phi_map, map, nside_factor)
   !AL: Added Oct 2007 
   !Temperature-only internally uses equi-cylindrical pix with bicubic interpolation
    use MPIstuff
    use Grid_Interpolation
    Type (HealpixInfo) :: H, H_res

    INTEGER(I4B), INTENT(IN) :: inlmax
    REAL(I4B), INTENT(IN), optional :: nside_factor
 
    integer nsmax
    real  :: nside_fac = 3.
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map
    REAL(SP), DIMENSION(:), pointer :: map2N, map2S
    
    Type (LensGradients) :: grad
    
    REAL(SP), DIMENSION(:,:), pointer :: high_res
     !Put two strips into same array so that at equator area is contiguous
     
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    INTEGER(I4B) :: n_phi_cyl,  cyl_start_ix, cyl_end_ix, cyl_start_edge, cyl_end_edge
    REAL(DP) :: dtheta_cyl ,dphi_cyl
    
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2,  corfac
    COMPLEX(DPC) :: factor, factor1
    COMPLEX(DPC) :: b_n, b_s

    CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2LENSEDMAPINTERPCYL'
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring,par_lm, nlmax

    COMPLEX(SPC) :: this_grad
    REAL(SP), DIMENSION(:), allocatable ::  ring, theta_vals, &
           phi_vals,theta_lensed_vals,phi_lensed_vals

    integer(I4B) a_ix, nalms, high_nside
    integer(I8B) ring_ix 
    REAL(DP) :: grad_len, sinc_grad_len, phi
    REAL(DP) :: cth0, sth0
    
    integer Ierror, InterpW
    real(I4B), dimension(:,:,:), allocatable :: Work
    integer(I4B) :: n_theta_cyl, cyl_th_end_ix 
    integer pole_edge, lmin
    integer status, topbottom
    COMPLEX(SPC), DIMENSION(:), pointer :: gradient_ring
#ifdef MPIPIX    
    double precision Initime
#endif  
    !=======================================================================

     nsmax = H%nside
     nlmax = inlmax
     if (present(nside_factor)) nside_fac = nside_factor

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if

    call SyncInts(nlmax)
    call SyncReals(nside_fac)
#endif
 
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms))
     if (H%MpiId==0) call Alm2PackAlm(alm,alm2,nlmax)
    
#ifdef MPIPIX 
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) then
      print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
      StartTime = geteTime()
     end if
#endif

!use equicylindrical high-res map section
    call cyl_interp_init(H, Grad, grad_phi_map, nside_fac,  high_nside, n_phi_cyl,dphi_cyl,dtheta_cyl, &
            cyl_start_ix, cyl_end_ix)

    call HealpixInitTrig(H_res,high_nside,nlmax, not_healpix = .true.)

    H_res%MpiId = 1
        
    cyl_start_ix = max(0,cyl_start_ix-interp_edge)
    cyl_end_ix = min(high_nside-1,cyl_end_ix+interp_edge )
    
    n_theta_cyl = cyl_end_ix - cyl_start_ix+1
    cyl_th_end_ix = cyl_end_ix+ n_theta_cyl
    
    if (cyl_start_ix ==0) then
      pole_edge = interp_edge !reflection about N/S pole for interpolation there
    else
      pole_edge = 0
    end if
    cyl_start_edge = cyl_start_ix - pole_edge
    cyl_end_edge = cyl_th_end_ix + pole_edge

    allocate(high_res(-interp_edge:n_phi_cyl-1+interp_edge, cyl_start_edge:cyl_end_edge), stat = status)
    if (status /= 0) call die_alloc(code,'high_res')
  
    ALLOCATE(ring(0:max(4*nsmax,n_phi_cyl)-1), stat = status)
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)

 
    do ith = cyl_start_ix, cyl_end_ix      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       !put pixels starting on phi=0, so kphi0=1 and first pixel centre at dtheta_cyl/2    
           
       cth = cos((ith+0.5_dp)*dtheta_cyl)
       sth = sin((ith+0.5_dp)*dtheta_cyl)
       kphi0 = 1
       
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       mmax_ring = get_mmax(nlmax,sth) 

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
         
          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW 
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac    !  actual lambda_mm      

          a_ix = a_ix + 1
          
          b_n = lam_lm * alm2(a_ix)
          b_s = b_n

          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          
          lmin = l_min_ylm(m,sth)
          
          do l = m+1, nlmax-1, 2
          !This is semi-optimized version where we do two at once
          
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             a_ix = a_ix + 1
         
             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

          !Second step unwind
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed

             if (l>lmin) then
              factor1 = lam_lm1m * alm2(a_ix)
              factor =  lam_lm * alm2(a_ix+1)
              b_n = b_n +  factor + factor1
              b_s = b_s +  factor - factor1
             end if
             lam_0 = lam_1 / a_rec
             lam_1 = lam_2

             a_ix = a_ix + 1
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo
          if (mod(nlmax - m,2)==1) then
          !do last one
             a_ix = a_ix + 1
             lam_lm = lam_2 * corfac 
             factor = lam_lm * alm2(a_ix)
             b_n = b_n + factor
             b_s = b_s - factor              
          end if
          

          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H_res,nlmax,b_north,n_phi_cyl,ring,kphi0,mmax_ring)   
       high_res(0:n_phi_cyl-1,ith) = ring(0:n_phi_cyl-1)
  
       call spinring_synthesis(H_res,nlmax, b_south, n_phi_cyl, ring, kphi0,mmax_ring)
       high_res(0:n_phi_cyl-1,cyl_th_end_ix - (ith-cyl_start_ix) ) = ring(0:n_phi_cyl-1)

    enddo    ! loop on cos(theta)

    deallocate(alm2)
    deallocate(H_res%trig)
    deallocate(H%recfac,stat= status)
    nullify(H%recfac)


#ifdef MPIPIX
    if(DebugMsgs>1) then
     print *,code//' Got high res ',H%MpiId, GeteTime()-StartTime
    end if 
    StartTime = Getetime()

     allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1), stat = status)
     if (status /= 0) call die_alloc(code,'map2N')
     allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1),stat = status)
     if (status /= 0) call die_alloc(code,'map2S')
#else
     map2N => map
     map2S => map
#endif

    allocate(phi_vals(-interp_edge:n_phi_cyl-1+interp_edge),theta_vals(cyl_start_edge:cyl_end_edge))
    do ring_ix = -interp_edge, n_phi_cyl-1 +interp_edge
      phi_vals(ring_ix) = (real(ring_ix,dp)+0.5_dp) * dphi_cyl
    end do
    if (pole_edge>0) then
      do ith = -pole_edge,-1
       !edge from over N/S pole for interpolation
        
        high_res(0:n_phi_cyl/2-1,ith) = high_res(n_phi_cyl-n_phi_cyl/2:n_phi_cyl-1 ,-ith-1) 
        high_res(n_phi_cyl/2:n_phi_cyl-1,ith) = high_res(0:n_phi_cyl-n_phi_cyl/2-1,-ith-1) 
        
        high_res(0:n_phi_cyl/2-1,cyl_th_end_ix-ith ) = high_res(n_phi_cyl-n_phi_cyl/2:n_phi_cyl-1,cyl_th_end_ix+ith+1)   
        high_res(n_phi_cyl/2:n_phi_cyl-1,cyl_th_end_ix-ith ) = high_res(0:n_phi_cyl-n_phi_cyl/2-1,cyl_th_end_ix+ith+1)   
    
   !    high_res(0:n_phi_cyl-1,ith) = high_res(n_phi_cyl-1:0:-1 ,-ith-1)   
   !    high_res(0:n_phi_cyl-1,cyl_th_end_ix-ith ) = high_res(n_phi_cyl-1:0:-1 ,cyl_th_end_ix+ith+1)   
      end do     
    end if
    !Patch in from both ends in phi so nice and smooth for interpolation 
    ! (can't be bothered to try to change interpolation for periodic boundary conditions)
    do ith = cyl_start_edge, cyl_end_edge
!     high_res(-interp_edge:-1,ith) = high_res(n_phi_cyl-1-interp_edge:n_phi_cyl-1-1,ith)
      high_res(-interp_edge:-1,ith) = high_res(n_phi_cyl-interp_edge:n_phi_cyl-1,ith)
      high_res(n_phi_cyl:n_phi_cyl+interp_edge-1,ith) = high_res(0:interp_edge-1,ith)
    end do
    
    do ith = cyl_start_edge,cyl_end_ix 
      theta_vals(ith) = (real(ith,dp)+0.5_dp) * dtheta_cyl
      theta_vals(cyl_th_end_ix - (ith-cyl_start_ix)) = pi -  theta_vals(ith)       
    end do
    

    allocate(phi_lensed_vals(0:nsmax*4-1),theta_lensed_vals(0:nsmax*4-1))

    allocate(Work(3,n_phi_cyl+interp_edge*2,cyl_start_edge:cyl_end_edge))
    InterpW=1 !Initialize interpolation
    
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

 
    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1

       if (ith < nsmax) then  ! polar cap (north)
          cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth0 = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
       endif
       
       do topbottom = 1,-1,-2
       
       if (topbottom==1) then
        gradient_ring => Grad%grad_phiN(H%istart_north(ith-1):)
       else
        gradient_ring => Grad%grad_phiS(H%istart_south(ith):)
       end if
       
       do ring_ix = 0, nph-1
    
        phi = ring_ix*2*pi/nph
!!        this_grad = Grad%grad_phiN(H%istart_north(ith-1)+ring_ix)
        this_grad = gradient_ring(ring_ix+1)
        grad_len = abs(this_grad)
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  topbottom*cos(grad_len) * cth0 - sinc_grad_len *sth0*real(this_grad) 
            sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
            phi = phi  +  asin(max(-1._dp,min(1._dp,aimag(this_grad)*sinc_grad_len/ sth )))
        else
         cth=topbottom*cth0
        endif

        if (kphi0==1) phi=phi + pi/nph
        if (phi >2._dp*pi) then
              phi = phi- 2._dp*pi  
           else if (phi< 0._dp) then
              phi = 2._dp*pi + phi
           end if 
           
        theta_lensed_vals(ring_ix) = dacos(cth)
        phi_lensed_vals(ring_ix) = phi

      end do

 !2D Cubic Interpolation using toms760 (smart version of bicubic) 
 !Most of the time in this stage is spent in here. ~50x more than calculating angles.
 !Mostly in the first call calculating grid of second derivatives
      call rgbi3p(Work,InterpW, n_phi_cyl+interp_edge*2, cyl_end_edge-cyl_start_edge+1, phi_vals, theta_vals, & 
          high_res, nph, phi_lensed_vals, theta_lensed_vals, ring, ierror)
 
! ITPLBV is toms474, and older routine. Much faster, but not as accute 
! and seems to be unstable to increasing the resolution
!       call ITPLBV(0, n_phi_cyl+interp_edge*2, cyl_end_edge-cyl_start_edge+1, phi_vals, theta_vals, & 
!          high_res, nph, phi_lensed_vals, theta_lensed_vals, ring) 
       
    
      InterpW = 2 !Don't need to re-compute next time as same grid
      if (ierror/=0) then
        write (*,*) code//': Interpolation error', Ierror
        stop
      end if  
      if (topBottom==1) then
          map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = ring(0:nph-1)
      else
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1) = ring(0:nph-1)
      end if    

     if (ith==2*nsmax) exit
     
     end do !loop on topBottom
 
     end do  !ith
 
    deallocate(ring)
    deallocate(Work)
  
    deallocate(theta_vals,phi_vals)
    deallocate(theta_lensed_vals,phi_lensed_vals)

    deallocate(high_res)
    call HealpixFree(H_res)

#ifdef MPIPIX
   
    deallocate(Grad%grad_phiN,Grad%grad_phiS)
   
    if(DebugMsgs>1) then
     print *,code//' Gather ',H%MpiId, GeteTime()-StartTime
    end if 
    StartTime = Getetime()
    call MPI_GATHERV(map2N(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),SP_MPI, &
       map(:),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2S(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),SP_MPI, &
       map(:),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2N,map2S)
     
#endif

  end subroutine scalalm2LensedmapInterpCyl


  subroutine alm2LensedmapInterpCyl(H,inlmax, alm_TEB, grad_phi_map, map_TQU, nside_factor)
    use MPIstuff
    use Grid_Interpolation

    Type (HealpixInfo) :: H, H_res

    INTEGER(I4B), INTENT(IN) :: inlmax
    REAL(I4B), INTENT(IN), optional :: nside_factor
 
    integer nsmax   
    real  :: nside_fac = 3.
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    REAL(SP), DIMENSION(:,:), pointer :: map2N, map2S

    Type (LensGradients) :: grad
    COMPLEX(SPC), DIMENSION(:), pointer :: gradient_ring
    REAL(SP), DIMENSION(:,:,:), pointer :: high_res
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB
    real(I4B), dimension(:,:,:,:), allocatable :: Work

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    INTEGER(I4B) :: n_phi_cyl,  cyl_start_ix, cyl_end_ix, cyl_start_edge, cyl_end_edge
    REAL(DP) :: dtheta_cyl ,dphi_cyl
    integer n_theta_cyl, cyl_th_end_ix,pole_edge

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fl,fm, f2m, fm2,  corfac, fm2fac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n, b_s
    COMPLEX(SPC) :: this_grad

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDMAPINTERPCYL'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact

    REAL(SP), DIMENSION(:), allocatable ::  ring, theta_vals, &
         phi_vals,theta_lensed_vals,phi_lensed_vals

    COMPLEX(SP), DIMENSION(:), allocatable ::  polrot

    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    real(dp), dimension(:), allocatable :: twocthlm1, lfac
   
    integer(I4B) a_ix, nalms, high_nside, interpW
    integer(I8B) ring_ix 
    COMPLEX(DPC) :: gammfac
    REAL(DP) :: grad_len, sinc_grad_len, phi, gamma 
    REAL(DP) :: Re, Im, cth0, sth0 
    integer lmin, polix, topbottom, ierror
#ifdef MPIPIX
    double precision Initime
    integer i
#endif
    !=======================================================================

    !     --- allocates space for arrays ---

     nsmax = H%nside
     nlmax = inlmax
     if (present(nside_factor)) nside_fac = nside_factor

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if

    call SyncInts(nlmax)
    call SyncReals(nside_fac)
#endif

     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(TEB(3,nalms))
     if (H%MpiId==0) then
        call TEB2PackTEB(alm_TEB,TEB,nlmax)
     end if
     
#ifdef MPIPIX
     call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) then
      print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
      StartTime = geteTime()
     end if
#endif


  !use equicylindrical high-res map section
    call cyl_interp_init(H, Grad, grad_phi_map, nside_fac,  high_nside, n_phi_cyl,dphi_cyl,dtheta_cyl, &
            cyl_start_ix, cyl_end_ix)

    call HealpixInitTrig(H_res,high_nside,nlmax, not_healpix = .true.)

    H_res%MpiId = 1
        
    cyl_start_ix = max(0,cyl_start_ix-interp_edge)
    cyl_end_ix = min(high_nside-1,cyl_end_ix+interp_edge )
    
    n_theta_cyl = cyl_end_ix - cyl_start_ix+1
    cyl_th_end_ix = cyl_end_ix+ n_theta_cyl
    
    if (cyl_start_ix ==0) then
      pole_edge = interp_edge !reflection about N/S pole for interpolation there
    else
      pole_edge = 0
    end if
    cyl_start_edge = cyl_start_ix - pole_edge
    cyl_end_edge = cyl_th_end_ix + pole_edge

    allocate(high_res(-interp_edge:n_phi_cyl-1+interp_edge, cyl_start_edge:cyl_end_edge,3), stat = status)
    if (status /= 0) call die_alloc(code,'high_res')
      
    ALLOCATE(ring(0:max(4*nsmax,n_phi_cyl)-1), stat = status)
    if (status /= 0) call die_alloc(code,'ring')

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')
    
    allocate(twocthlm1(0:nlmax))
    allocate(lfac(nlmax))

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')


    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   lam_fact = lam_fact * 2 !HealPix polarization def

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
       normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(high_nside)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(high_nside))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(high_nside) )


     do ith = cyl_start_ix, cyl_end_ix      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       !put pixels starting on phi=0, so kphi0=1 and first pixel centre at dtheta_cyl/2    
           
           
       cth = cos((ith+0.5_dp)*dtheta_cyl)
       sth = sin((ith+0.5_dp)*dtheta_cyl)
       kphi0 = 1
       
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       mmax_ring = get_mmax(nlmax,sth) 

       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2

        do l=1,nlmax
        twocthlm1(l) = cth*real(2*(l-1),dp) !! 2(l-1)cos(theta)
        lfac(l) = -real(2*l,dp)*one_on_s2 -real(l,dp)*real(l-1,dp)
       end do

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          fm2fac= 2._dp* fm2 * one_on_s2
          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac     !  actual lambda_mm      

          a_ix = a_ix + 1

          b_n = lam_lm * TEB(1,a_ix)
          b_s = b_n

          !l=m special case
          if (m >=2) then
              lambda_w = - 2.0_dp *(normal_l(m) * lam_lm * (fm - fm2) ) * ( one_on_s2 - 0.5_dp )
              lambda_x =  ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              b_n_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_s_Q =  par_lm*(lambda_w * TEB(2,a_ix) - zi_lam_x * TEB(3,a_ix))

              b_n_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
              b_s_U = par_lm*(lambda_w * TEB(3,a_ix) + zi_lam_x * TEB(2,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec

          lmin = max(2,l_min_ylm(m,sth))

          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             a_ix = a_ix + 1

             factor = lam_lm * TEB(1,a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

             if (l>=lmin .and. corfac /= 0) then

                 a_w =   fm2fac  + lfac(l) !2* (fm2 - fl) * one_on_s2 - fl*(fl - 1._dp)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  twocthlm1(l) * lam_lm      !2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U
          b_north(m) = b_n
          b_south(m) = b_s

       enddo
!North
       call spinring_synthesis(H_res,nlmax,b_north,n_phi_cyl,ring,kphi0,mmax_ring)   
       high_res(0:n_phi_cyl-1,ith,1) = ring(0:n_phi_cyl-1)
       call spinring_synthesis(H_res,nlmax, b_north_Q, n_phi_cyl, ring, kphi0,mmax_ring)
       high_res(0:n_phi_cyl-1,ith,2) = ring(0:n_phi_cyl-1)
       call spinring_synthesis(H_res,nlmax, b_north_U, n_phi_cyl, ring, kphi0,mmax_ring)
       high_res(0:n_phi_cyl-1,ith,3) = ring(0:n_phi_cyl-1)
!South  
       call spinring_synthesis(H_res,nlmax, b_south, n_phi_cyl, ring, kphi0,mmax_ring)
       high_res(0:n_phi_cyl-1,cyl_th_end_ix - (ith-cyl_start_ix) ,1) = ring(0:n_phi_cyl-1)
       call spinring_synthesis(H_res,nlmax, b_south_Q, n_phi_cyl, ring, kphi0,mmax_ring)
       high_res(0:n_phi_cyl-1,cyl_th_end_ix - (ith-cyl_start_ix) ,2 ) = ring(0:n_phi_cyl-1)
       call spinring_synthesis(H_res,nlmax, b_south_U, n_phi_cyl, ring, kphi0,mmax_ring)
       high_res(0:n_phi_cyl-1,cyl_th_end_ix - (ith-cyl_start_ix) ,3) = ring(0:n_phi_cyl-1)

    enddo    ! loop on cos(theta)


    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l,twocthlm1,lfac)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)
    deallocate(TEB)
    deallocate(H_res%trig)

    deallocate(H%recfac,stat= status)
    nullify(H%recfac)

#ifdef MPIPIX
    if(DebugMsgs>1) then
     print *,code//' Got high res ',H%MpiId, GeteTime()-StartTime
    end if 
    StartTime = Getetime()

     allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,3), stat = status)
     if (status /= 0) call die_alloc(code,'map2N')
     allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,3),stat = status)
     if (status /= 0) call die_alloc(code,'map2S')
#else
     map2N => map_TQU
     map2S => map_TQU
#endif

    allocate(phi_vals(-interp_edge:n_phi_cyl-1+interp_edge),theta_vals(cyl_start_edge:cyl_end_edge))
    do ring_ix = -interp_edge, n_phi_cyl-1 +interp_edge
      phi_vals(ring_ix) = (real(ring_ix,dp)+0.5_dp) * dphi_cyl
    end do
    if (pole_edge>0) then
      do ith = -pole_edge,-1
       !edge from over N/S pole for interpolation
       high_res(0:n_phi_cyl/2-1,ith,:) = high_res(n_phi_cyl-n_phi_cyl/2:n_phi_cyl-1 ,-ith-1,:) 
       high_res(n_phi_cyl/2:n_phi_cyl-1,ith,:) = high_res(0:n_phi_cyl-n_phi_cyl/2-1,-ith-1,:) 
        
       high_res(0:n_phi_cyl/2-1,cyl_th_end_ix-ith,:) = high_res(n_phi_cyl-n_phi_cyl/2:n_phi_cyl-1,cyl_th_end_ix+ith+1,:)   
       high_res(n_phi_cyl/2:n_phi_cyl-1,cyl_th_end_ix-ith,:) = high_res(0:n_phi_cyl-n_phi_cyl/2-1,cyl_th_end_ix+ith+1,:)   
  
      end do     
    end if
    !Patch in from both ends in phi so nice and smooth for interpolation 
    ! (can't be bothered to try to change interpolation for periodic boundary conditions)
    do ith = cyl_start_edge, cyl_end_edge
     high_res(-interp_edge:-1,ith,:) = high_res(n_phi_cyl-interp_edge:n_phi_cyl-1,ith,:)
     high_res(n_phi_cyl:n_phi_cyl+interp_edge-1,ith,:) = high_res(0:interp_edge-1,ith,:)
    end do
    
    do ith = cyl_start_edge,cyl_end_ix 
      theta_vals(ith) = (real(ith,dp)+0.5_dp) * dtheta_cyl
      theta_vals(cyl_th_end_ix - (ith-cyl_start_ix)) = pi -  theta_vals(ith)       
    end do
    

    ALLOCATE(polrot(0:4*nsmax-1), stat = status)

    allocate(phi_lensed_vals(0:nsmax*4-1),theta_lensed_vals(0:nsmax*4-1))

    allocate(Work(3,n_phi_cyl+interp_edge*2,cyl_start_edge:cyl_end_edge,3))
    InterpW=1 !Initialize interpolation
    

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1

       if (ith < nsmax) then  ! polar cap (north)
          cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth0 = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
       endif

       do topbottom = 1,-1,-2 !+1 or -1 for north and south

       if (topbottom==1) then
        gradient_ring => Grad%grad_phiN(H%istart_north(ith-1):)
       else
        gradient_ring => Grad%grad_phiS(H%istart_south(ith):)
       end if
       
       do ring_ix = 0, nph-1

        phi = ring_ix*2*pi/nph
        this_grad = gradient_ring(ring_ix+1)
        grad_len = abs(this_grad)
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  topbottom*cos(grad_len) * cth0 - sinc_grad_len *sth0*real(this_grad) 
            sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
            phi = phi  + asin(max(-1._dp,min(1._dp,aimag(this_grad)*sinc_grad_len/ sth )))
        else
         cth=topbottom*cth0
         sth=sth0
        endif

        if (kphi0==1) phi=phi + pi/nph
        if (phi >2._dp*pi) then
              phi = phi- 2._dp*pi  
           else if (phi< 0._dp) then
              phi = 2._dp*pi + phi
           end if

        theta_lensed_vals(ring_ix) = dacos(cth)
        phi_lensed_vals(ring_ix) = phi

         if (grad_len > 1e-20_dp) then
             Re = real(this_grad)
             Im = aimag(this_grad)
             gamma = topbottom*grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             polrot(ring_ix)=  cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
         else 
           polrot(ring_ix)=1
         end if

      end do

 !2D Cubic Interpolate 
 !Most of the time in this stage is spent in here. ~50x more than calculating angles
 !Should probably try something faster
     do polix =1,3
      call rgbi3p(Work(:,:,:,polix),InterpW, n_phi_cyl+interp_edge*2, cyl_end_edge-cyl_start_edge+1, phi_vals, theta_vals, & 
          high_res(:,:,polix), nph, phi_lensed_vals, theta_lensed_vals, ring, ierror)
    
      if (ierror/=0) then
        write (*,*) code // ': Interpolation error', Ierror
        stop
      end if  
      if (topBottom==1) then
          map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,polix) = ring(0:nph-1)
      else
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1,polix) = ring(0:nph-1)
      end if    
     end do
     InterpW = 2 !Don't need to re-compute next time as same grid

     !Put in rotation factor
      do ring_ix = 0, nph-1
        if (topbottom==1) then    
         gammfac = polrot(ring_ix)*cmplx(map2N(H%istart_north(ith-1)+ring_ix,2),map2N(H%istart_north(ith-1)+ring_ix,3))           
         map2N(H%istart_north(ith-1)+ring_ix,2) = real(gammfac)
         map2N(H%istart_north(ith-1)+ring_ix,3) = aimag(gammfac)
        else
         gammfac = polrot(ring_ix)*cmplx(map2S(H%istart_south(ith)+ring_ix,2),map2S(H%istart_south(ith)+ring_ix,3))           
         map2S(H%istart_south(ith)+ring_ix,2) = real(gammfac)
         map2S(H%istart_south(ith)+ring_ix,3) = aimag(gammfac)
       endif
      end do

     if (ith==2*nsmax) exit

     end do !loop on topBottom

    end do !ith


    deallocate(ring)
    deallocate(Work)
  
    deallocate(theta_vals,phi_vals)
    deallocate(theta_lensed_vals,phi_lensed_vals)

    deallocate(high_res)
    call HealpixFree(H_res)

    deallocate(polrot)


#ifdef MPIPIX

    deallocate(Grad%grad_phiN,Grad%grad_phiS)
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiId, GeteTime()-StartTime
    StartTime = Getetime()
    do i=1,3
    call MPI_GATHERV(map2N(H%North_Start(H%MpiId),i),H%North_Size(H%MpiId),SP_MPI, &
       map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2S(H%South_Start(H%MpiId),i),H%South_Size(H%MpiId),SP_MPI, &
       map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end do
    if(DebugMsgs>1) print *,code //' Wait at gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2N,map2S) 
#endif

  end subroutine alm2LensedmapInterpCyl



  subroutine alm2LensedQuadContrib(H, inlmax, alm, grad_phi_map, map_T)
 !Get lensed map using fourth order series expansion
 !*** Note this does not give a map with an accurate lensed C_l at l>~ 1200***
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H
    INTEGER(I4B), INTENT(IN) :: inlmax 
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:)  :: alm
    REAL(SPC), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map_T
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SPC), DIMENSION(:), pointer :: map2
    COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0, nlmax

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: a_w, b_w, a_x, fac
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n, b_s, b_n4,b_s4     
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n_Q3, b_s_Q3, b_n_U3, b_s_U3
    COMPLEX(DPC) :: b_n_Q4, b_s_Q4, b_n_U4, b_s_U4
    COMPLEX(DPC) :: b_n_Q33, b_s_Q33, b_n_U33, b_s_U33
    COMPLEX(DPC) :: b_n_Q44, b_s_Q44, b_n_U44, b_s_U44
    COMPLEX(DPC) :: b_n_Q2, b_s_Q2, b_n_U2, b_s_U2

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDQUADCONTRIB'
    COMPLEX(DPC) ::  b(0:H%lmax,2), b4(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q(0:H%lmax,2), b_U(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q2(0:H%lmax,2), b_U2(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q3(0:H%lmax,2), b_U3(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q33(0:H%lmax,2), b_U33(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q44(0:H%lmax,2), b_U44(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q4(0:H%lmax,2), b_U4(0:H%lmax,2)

    REAL(DP) :: LastX(4),LastW(4),W(4),X(4)
    
    INTEGER(I4B) :: status,par_lm, a_ix

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:),   ALLOCATABLE :: ringR, ringI
    integer NS, mmax_ring, nalms, ix
#ifdef MPIPIX    
    double precision Initime
#endif      !=======================================================================

     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call SyncInts(nlmax)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms),stat = status )
     if (status /= 0) call die_alloc(code,'alm2')
     if (H%MpiId==0) then
       call Alm2PackAlm(alm,alm2,nlmax)
       grad_phi => grad_phi_map
     else
       allocate(grad_phi(0:12*H%nside**2-1)) 
     end if
#ifdef MPIPIX
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code//' Got alm ',H%MpiId, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)
     if (status /= 0) call die_alloc(code,'map2')    
#else
     map2 => map_T
#endif

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(ringR(0:4*nsmax-1),ringI(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    call HealpixInitRecfac(H,nlmax)
    call GetLamfact(lam_fact, nlmax)
!Don't put in 2 factor for spin 2 here, but below   
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )
    !     --------------------------------------------

    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)   ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       
       mmax_ring = get_mmax(nlmax,sth) 

       lam_mm = sq4pi_inv ! lamda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2
          !           ---------- l = m ----------
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW

          lam_lm = corfac     !  actual lambda_mm      
          
          a_ix = a_ix + 1

          par_lm = 1  

!Del^2 T term

          b_n = -(fm2+fm)*lam_lm * alm2(a_ix)
          b_s = b_n

!4th order
          fac = (2-3*(fm2+fm))
          b_n4 = factor * b_n
          b_s4 = b_n4

          par_lm = -1
!Grad T term
          if (m >=1) then
!normal_l cancels with gradient, sign from gradient 
              X(1) =  lam_lm * fm / sth 
              W(1) =  -X(1) * cth

              b_n_Q =  -W(1) * alm2(a_ix)
              b_s_Q =  par_lm * b_n_Q

              b_n_U = cmplx(0._dp,X(1)) * alm2(a_ix)
              b_s_U = -par_lm * b_n_U

              b_n_Q3 =  -fac * W(1) * alm2(a_ix)
              b_s_Q3 =  par_lm * b_n_Q

              b_n_U3 = cmplx(0._dp,X(1)*fac) * alm2(a_ix)
              b_s_U3 = -par_lm * b_n_U

          else
             W(1) = 0
             X(1) = 0

             b_n_Q=0; b_s_Q=0; b_n_U=0; b_s_U=0
             b_n_Q3=0; b_s_Q3=0; b_n_U3=0; b_s_U3=0
          end if

!eth eth T term
          par_lm = 1
          if (m >=2) then
              W(2) = -( lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
              X(2) =  ( lam_lm * (fm - fm2) ) *   2.0_dp * c_on_s2
              b_n_Q2 =  W(2)* alm2(a_ix) 
              b_s_Q2 =  par_lm* b_n_Q2

              b_n_U2 =  cmplx(0.0_dp, -X(2)) * alm2(a_ix)
              b_s_U2 =  -par_lm*b_n_U2  
 
 !4th order spin 2 term
              b_n_Q4 =  (8 - 4*(fm+fm2)) * W(2)* alm2(a_ix) 
              b_s_Q4 =  par_lm* b_n_Q4

              b_n_U4 =  cmplx(0.0_dp, -(8 - 4*(fm+fm2))*X(2)) * alm2(a_ix)
              b_s_U4 =  -par_lm*b_n_U4  

          else
             W(2)=0;X(2)=0
             b_n_Q2=0;b_s_Q2=0;b_n_U2=0;b_s_U2=0
             b_n_Q4=0;b_s_Q4=0;b_n_U4=0;b_s_U4=0
          end if

!Irreducible spin 3 term
          par_lm = -1
          if (m >=3) then
              X(3) =  lam_lm / sth * fm*(fm-1)*(fm-2) 
              W(3) = X(3) * cth * ( 1 - 4*one_on_s2)
              X(3) = X(3) * (4*one_on_s2 - 3)             

              b_n_Q33 =  W(3) * alm2(a_ix) 
              b_s_Q33 =  par_lm* b_n_Q2

              b_n_U33 =  cmplx(0.0_dp, -X(3)) * alm2(a_ix)
              b_s_U33 =  -par_lm*b_n_U2  

          else
             W(3)=0;X(3)=0
             b_n_Q33=0;b_s_Q33=0;b_n_U33=0;b_s_U33=0
          end if

          par_lm = 1
  
           if (m>=4 ) then
   
             W(4) = ((l-3)*(X(3) - cth*W(3))) /(sth)
             X(4) = ((l-3)*(W(3) - cth*X(3))) /(sth)
             factor_1 =  W(4) * alm2(a_ix)
             b_n_Q44 =  factor_1
             b_s_Q44 = par_lm * factor_1 

             factor_2 =  cmplx(0.0_dp, -X(4)) * alm2(a_ix) 
             b_n_U44 =  factor_2
             b_s_U44 =  -par_lm * factor_2
           else 
             b_n_Q44=0;b_s_Q44=0;b_n_U44=0;b_s_U44=0
           end if


!Keep par_lm correct for spin zero
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m
             LastX = X
             LastW = W 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl

             a_ix = a_ix + 1

!Del^2 T
             factor = -(fl2+fl)*lam_lm * alm2(a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

!4th order
             fac = (2-3*(fl2+fl))
             factor = factor * fac 
             b_n4 = b_n4 + factor
             b_s4 = b_s4 + par_lm*factor

!grad T
             par_lm = -par_lm
             a_w = 1 / sth
             X(1) = a_w * fm  * lam_lm
             W(1) = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 

             factor_1 =  -W(1) * alm2(a_ix)
             b_n_Q = b_n_Q +          factor_1 
             b_s_Q = b_s_Q + par_lm * factor_1 ! X has a diff. parity

             factor_2 =   cmplx(0._dp,-X(1)) * alm2(a_ix) 
             b_n_U = b_n_U - factor_2
             b_s_U = b_s_U + par_lm * factor_2

!grad [(3 grad^2 + 2)T]
             factor_1 = - fac * W(1) * alm2(a_ix)
             b_n_Q3 = b_n_Q3 +          factor_1 
             b_s_Q3 = b_s_Q3 + par_lm * factor_1 ! X has a diff. parity

             factor_2 =   cmplx(0._dp,-fac * X(1)) * alm2(a_ix) 
             b_n_U3 = b_n_U3 - factor_2
             b_s_U3 = b_s_U3 + par_lm * factor_2


!eth eth T
             par_lm = -par_lm
             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
    !put in 2 factor here
                 b_w =  2* c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 W(2) =  ( a_w * lam_lm + b_w * lam_lm1m ) 
    !and here
                 X(2) =  fm_on_s2 * ( 2* lam_fact(a_ix) * lam_lm1m - a_x)

                 factor_1 =  W(2)* alm2(a_ix)
                 b_n_Q2 = b_n_Q2 +  factor_1
                 b_s_Q2 = b_s_Q2 + par_lm * factor_1 

                 factor_2 =  cmplx(0.0_dp, X(2)) * alm2(a_ix) 
                 b_n_U2 = b_n_U2  - factor_2
                 b_s_U2 = b_s_U2 + par_lm * factor_2
    
    !4th order spin 2 term

                 factor_1 = (8-4*(fl2+fl))* W(2)* alm2(a_ix)
                 b_n_Q4 = b_n_Q4 +  factor_1
                 b_s_Q4 = b_s_Q4 + par_lm * factor_1 

                 factor_2 =  cmplx(0.0_dp, (8-4*(fl2+fl))* X(2)) * alm2(a_ix) 
                 b_n_U4 = b_n_U4  - factor_2
                 b_s_U4 = b_s_U4 + par_lm * factor_2

              else
               X(2)=0;W(2)=0;
              end if
!Irreducible spin 3 term
             par_lm = -par_lm
             if (l>=3 .and. corfac /= 0) then
   
           W(3) = ((l+2)*lam_fact(a_ix)*LastW(2) + (l-2)*(m*X(2) - cth*l*W(2))) /(l*sth)
           X(3) = ((l+2)*lam_fact(a_ix)*LastX(2) + (l-2)*(m*W(2) - cth*l*X(2))) /(l*sth)
!              a_w = 1._dp /sth
!                 b_w = (l-1)*(l-2)
!                 lambda_x = a_w*fm*( (one_on_s2*(4*fm2-(12*l-8)) - 3 * b_w)*lam_lm + &
!                         12*lam_fact(a_ix) * c_on_s2  * lam_lm1m)  
 !                lambda_w = a_w*( (l*b_w - one_on_s2*(8*l+fm2*(4*l-12))) * cth * lam_lm  &
 !                       - lam_fact(a_ix)*( fl2 + fl +6 - (8+4*fm2)*one_on_s2) * lam_lm1m)
                 factor_1 =  W(3) * alm2(a_ix)
                 b_n_Q33 = b_n_Q33 +  factor_1
                 b_s_Q33 = b_s_Q33 + par_lm * factor_1 

                 factor_2 =  cmplx(0.0_dp, X(3)) * alm2(a_ix) 
                 b_n_U33 = b_n_U33  - factor_2
                 b_s_U33 = b_s_U33 + par_lm * factor_2
             else
              W(3)=0;X(3)=0
             end if
             par_lm = -par_lm

!Irreducible spin 4 term
           if (l>=4 .and. corfac /= 0) then
   
             W(4) = ((l+3)*lam_fact(a_ix)*LastW(3) + (l-3)*(m*X(3) - cth*l*W(3))) /(l*sth)
             X(4) = ((l+3)*lam_fact(a_ix)*LastX(3) + (l-3)*(m*W(3) - cth*l*X(3))) /(l*sth)
             factor_1 =  W(4) * alm2(a_ix)
             b_n_Q44 = b_n_Q44 +  factor_1
             b_s_Q44 = b_s_Q44 + par_lm * factor_1 

             factor_2 =  cmplx(0.0_dp, X(4)) * alm2(a_ix) 
             b_n_U44 = b_n_U44  - factor_2
             b_s_U44 = b_s_U44 + par_lm * factor_2
     
           end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo

          b(m,1) = b_n; b(m,2) = b_s
          b4(m,1) = b_n4; b4(m,2) = b_s4

          b_Q(m,1) = b_n_Q; b_Q(m,2) = b_s_Q; b_U(m,1) = b_n_U; b_U(m,2) = b_s_U
          b_Q2(m,1) = b_n_Q2; b_Q2(m,2) = b_s_Q2; b_U2(m,1) = b_n_U2; b_U2(m,2) = b_s_U2
          b_Q3(m,1) = b_n_Q3; b_Q3(m,2) = b_s_Q3; b_U3(m,1) = b_n_U3; b_U3(m,2) = b_s_U3
          b_Q4(m,1) = b_n_Q4; b_Q4(m,2) = b_s_Q4; b_U4(m,1) = b_n_U4; b_U4(m,2) = b_s_U4
          b_Q33(m,1) = b_n_Q33; b_Q33(m,2) = b_s_Q33; b_U33(m,1) = b_n_U33; b_U33(m,2) = b_s_U33
          b_Q44(m,1) = b_n_Q44; b_Q44(m,2) = b_s_Q44; b_U44(m,1) = b_n_U44; b_U44(m,2) = b_s_U44

       enddo

       ix = H%istart_north(ith-1)
       do NS = 1,2  
!Fourth order
   !spin zero
           call spinring_synthesis(H,nlmax, b4(:,NS), nph, ringR, kphi0,mmax_ring)
           map2(ix:ix+nph-1) =  &
           (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2)**2*RingR(0:nph-1)
            !factor 1/(4!) goes in below
   !irreducible spin 4
           call spinring_synthesis(H,nlmax, b_Q44(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U44(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + real( &
             grad_phi(ix:ix+nph-1)**4*cmplx(RingR(0:nph-1),-RingI(0:nph-1)) )
     
   !spin 2 cross
           call spinring_synthesis(H,nlmax, b_Q4(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U4(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = (map2(ix:ix+nph-1) + &
             (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2) * &
             real( grad_phi(ix:ix+nph-1)**2*cmplx(RingR(0:nph-1),-RingI(0:nph-1)) ))/8


!Cubic reducible term
           call spinring_synthesis(H,nlmax, b_Q3(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U3(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + ((real(grad_phi(ix:ix+nph-1))*RingR(0:nph-1) + &
                              aimag(grad_phi(ix:ix+nph-1))*RingI(0:nph-1)) * &
                       (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2))
!Cubic irreducible term
           call spinring_synthesis(H,nlmax, b_Q33(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U33(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = (map2(ix:ix+nph-1) + real( &
             -grad_phi(ix:ix+nph-1)**3*cmplx(RingR(0:nph-1),-RingI(0:nph-1)) )) / (3*2)
                   !factor of 1/4 goes in with quadratic addition


    !Quadratic
      ! Re(eth eth T ethb phi ethb phi)
           call spinring_synthesis(H,nlmax, b_Q2(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U2(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + &
                   (real(grad_phi(ix:ix+nph-1))+ aimag(grad_phi(ix:ix+nph-1)))* &
                   (real(grad_phi(ix:ix+nph-1))- aimag(grad_phi(ix:ix+nph-1))) * RingR(0:nph-1) &
              + 2*RingI(0:nph-1) * real(grad_phi(ix:ix+nph-1))*aimag(grad_phi(ix:ix+nph-1)) 
      !grad^2 T |Grad phi|^2
           call spinring_synthesis(H,nlmax, b(:,NS), nph, ringR, kphi0,mmax_ring)
           map2(ix:ix+nph-1) =0.25_dp * ( map2(ix:ix+nph-1) + &
               (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2)*RingR(0:nph-1))

    !Linear
     !grad T dot grad phi
           call spinring_synthesis(H,nlmax, b_Q(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + real(grad_phi(ix:ix+nph-1))*RingR(0:nph-1) + &
                              aimag(grad_phi(ix:ix+nph-1))*RingI(0:nph-1)

           if (ith  >=  2*nsmax) exit
           ix = H%istart_south(ith)
      end do

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(ringR,ringI)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiId
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiId)),H%North_Size(H%MpiId),SP_MPI, &
       map_T,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiId)),H%South_Size(H%MpiId),SP_MPI, &
       map_T,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if (H%MpiId/=0) deallocate(grad_phi)
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine alm2LensedQuadContrib



  subroutine map2polalm(h, inlmax, map_TQU,alm_TEB, cos_theta_cut)
    use mpistuff
    type (healpixinfo) :: h

    integer(i4b)  :: nsmax
    integer(i4b), intent(in) :: inlmax
    complex(spc), intent(out),  dimension(:,:,:) :: alm_TEB
    real(sp), intent(in), dimension(0:12*h%nside**2-1,3), target :: map_TQU
    complex(spc), dimension(:,:), allocatable :: TEB
    real(sp), dimension(:,:), pointer :: map2S,map2N

    real(dp),     intent(in) :: cos_theta_cut

    integer :: nlmax !base integer type for mpi compat
    integer(i4b) :: l, m, ith, scalem, scalel        ! alm related
    integer(i4b) :: nph, kphi0 ! map related

    real(dp) :: omega_pix
    real(dp) :: cth, sth, dth1, dth2, dst1
    real(dp) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    real(dp) :: fl, fm, f2m, fm2, fm2fac,corfac
    real(dp) :: c_on_s2, fm_on_s2, one_on_s2
    real(dp) :: lambda_w, lambda_x, a_w, b_w, a_x
  
    character(len=*), parameter :: code = 'MAP2POLALM'
    complex(dpc), dimension(:), allocatable :: phas_nq, phas_nu
    complex(dpc), dimension(:), allocatable :: phas_sq, phas_su
    complex(dpc), dimension(:), allocatable :: phas_n
    complex(dpc), dimension(:), allocatable :: phas_s
    complex(dpc) :: phas_Qp,phas_Qm,phas_UM,phas_Up,phas_p,phas_m
    complex(dpc) :: Iphas_Qp,Iphas_Qm,Iphas_UM,Iphas_Up
    

    integer(i4b) mmax_ring, status, par_lm, a_ix
    integer nalms, lmin
    real(dp), dimension(:), allocatable :: lam_fact
    real(dp), dimension(:),   allocatable :: ring
    real(dp), dimension(:),   allocatable :: normal_l
    real(dp), dimension(:), allocatable :: twocthlm1, lfac
    
    logical   :: keep_it
#ifdef MPIPIX
    double precision initime
    integer i
#endif
    !=======================================================================

    nsmax = h%nside
    nlmax = inlmax
     
#ifdef MPIPIX
     if (cos_theta_cut/=-1) stop 'cos_theta_cut /= -1'
     if (h%MpiId==0) then 
      if(debugmsgs>0) print *,code //': sending to farm '
      call sendmessages(h,code)
      map2S => map_TQU
      map2N => map_TQU
    else
       allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,3),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
       allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,3),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
    end if

    starttime = getetime()    
    initime = starttime
    call SyncInts(nlmax)
    do i=1,3
    call mpi_scatterv(map_TQU(:,i),h%north_size, h%north_start, &
       sp_mpi, map2N(h%north_start(h%MpiId),i),h%north_size(h%MpiId),sp_mpi, 0 ,mpi_comm_world, ierr)
    call mpi_scatterv(map_TQU(:,i),h%south_size, h%south_start, &
       sp_mpi, map2S(h%south_start(h%MpiId),i),h%south_size(h%MpiId),sp_mpi, 0 ,mpi_comm_world, ierr)
    end do
    if(debugmsgs>1) print *,code //' scattered ',h%MpiId, getetime() - starttime
#else
    map2N => map_TQU
    map2S => map_TQU
#endif

    nalms = ((nlmax+1)*(nlmax+2))/2   

    allocate(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    allocate(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')
    allocate(twocthlm1(0:nlmax))
    allocate(lfac(nlmax))
    allocate(phas_n(0:nlmax), phas_nq(0:nlmax),&
         &   phas_nu(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    allocate(phas_s(0:nlmax), phas_sq(0:nlmax),&
         &   phas_su(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')

    allocate(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

   call healpixinitrecfac(h,nlmax)
   call getlamfact(lam_fact, nlmax)
   lam_fact = lam_fact*2

   allocate(TEB(3,nalms),stat = status)
   if (status /= 0) call die_alloc(code,'TEB')
   TEB = 0    
       
   omega_pix = pi / (3 * nsmax * real(nsmax,dp))

    normal_l = 0.0_dp
    do l = 2, nlmax
        fl = dble(l)
        normal_l(l) = eb_sign * sqrt( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo

    dth1 = 1.0_dp / (3.0_dp*dble(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*dble(nsmax))
    dst1 = 1.0_dp / (sqrt(6.0_dp) * dble(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    
    do ith = h%ith_start(h%MpiId), h%ith_end(h%MpiId)

       phas_nq=0; phas_sq=0;phas_nu=0;phas_su=0; phas_n=0; phas_s=0

       if (ith  <=  nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - dble(ith)**2 * dth1
          sth = sin( 2.0_dp * asin( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = mod(ith+1-nsmax,2)
          cth = dble(2*nsmax-ith) * dth2
          sth = dsqrt((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       do l=1,nlmax
        twocthlm1(l) = cth*real(2*(l-1),dp) !! 2(l-1)cos(theta)
        lfac(l) = -real(2*l,dp)*one_on_s2 -real(l,dp)*real(l-1,dp)
       end do

       mmax_ring = get_mmax(nlmax,sth) 

       keep_it = (abs(cth) > cos_theta_cut) ! part of the sky out of the symmetric cut
       if (keep_it) then
          ring(0:nph-1) = map2N(h%istart_north(ith-1):h%istart_north(ith-1)+nph-1,1) * h%w8ring_TQU(ith,1)
          call spinring_analysis(h,nlmax, ring, nph, phas_n, kphi0, mmax_ring)
          ring(0:nph-1) = map2N(h%istart_north(ith-1):h%istart_north(ith-1)+nph-1,2) * h%w8ring_TQU(ith,2)
          call spinring_analysis(h,nlmax, ring, nph, phas_nq, kphi0, mmax_ring)
          ring(0:nph-1) = map2N(h%istart_north(ith-1):h%istart_north(ith-1)+nph-1,3) * h%w8ring_TQU(ith,3) 
          call spinring_analysis(h,nlmax, ring, nph, phas_nu, kphi0, mmax_ring)
       endif

       if (ith  <  2*nsmax .and. keep_it) then
          ring(0:nph-1) = map2S(h%istart_south(ith):h%istart_south(ith)+nph-1,1) * h%w8ring_TQU(ith,1)
          call spinring_analysis(h,nlmax, ring, nph, phas_s, kphi0, mmax_ring)
          ring(0:nph-1) = map2S(h%istart_south(ith):h%istart_south(ith)+nph-1,2) * h%w8ring_TQU(ith,2)
          call spinring_analysis(h,nlmax, ring, nph, phas_sq, kphi0, mmax_ring)
          ring(0:nph-1) = map2S(h%istart_south(ith):h%istart_south(ith)+nph-1,3) * h%w8ring_TQU(ith,3)
          call spinring_analysis(h,nlmax, ring, nph, phas_su, kphi0, mmax_ring)
       endif

       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       if (keep_it) then ! avoid un-necessary calculations (eh, 09-2001)
          lam_mm = sq4pi_inv * omega_pix 
          scalem=1
          a_ix = 0
          do m = 0, mmax_ring
             fm  = dble(m)
             f2m = 2.0_dp * fm
             fm2 = fm * fm
             
             fm2fac= 2._dp* fm2 * one_on_s2
             
             fm_on_s2 = fm * one_on_s2

             !           ---------- l = m ----------
             par_lm = 1   ! = (-1)^(l+m+s)
             if (m  >=  1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm) < unflow) then
                lam_mm=lam_mm*ovflow
                scalem=scalem-1
             endif

             a_ix = a_ix+1

             corfac = scalefactor(scalem)*lam_mm/ovflow
             lam_lm = corfac
          
             lmin = max(2,l_min_ylm(m, sth))
             phas_Qp=phas_nQ(m) + phas_sQ(m)
             phas_Qm=phas_nQ(m) - phas_sQ(m)
             phas_Up=phas_nU(m) + phas_sU(m)
             phas_Um=phas_nU(m) - phas_sU(m)
             Iphas_Qp = (0.0_dp, 1._dp)*phas_Qp
             Iphas_Qm = (0.0_dp, 1._dp)*phas_Qm
             Iphas_Up = (0.0_dp, 1._dp)*phas_Up
             Iphas_Um = (0.0_dp, 1._dp)*phas_Um

             phas_p= phas_n(m) + phas_s(m)
             phas_m= phas_n(m) - phas_s(m)


             TEB(1,a_ix) = TEB(1,a_ix) + lam_lm * phas_p
             if (m >= lmin) then
            
              lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) *( 2.0_dp * one_on_s2 - 1.0_dp )
              lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
                 TEB(2,a_ix) = TEB(2,a_ix) &
                  &                 + lambda_w * phas_Qp + lambda_x*Iphas_Um

                 TEB(3,a_ix) = TEB(3,a_ix) &
                  &                 + lambda_w * phas_Up - lambda_x*Iphas_Qm
            
             end if

             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = h%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec
         
             do l = m+1, nlmax
                par_lm = - par_lm  ! = (-1)^(l+m)
                lam_lm1m=lam_lm ! actual lambda_l-1,m (useful for polarisation)
                lam_lm   = lam_2*corfac ! actual lambda_lm (ovflow factors removed)
                
                a_ix = a_ix + 1

                TEB(1,a_ix) = TEB(1,a_ix) + lam_lm * (phas_n(m) + par_lm*phas_s(m))

             if (l>=lmin .and. corfac /= 0) then
                 !corfac=0 guarantees lam(l-1) is also v close to zero

!                a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)

                 a_w =   fm2fac  + lfac(l)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  twocthlm1(l) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)

       
                if (par_lm==1) then

                 TEB(2,a_ix) = TEB(2,a_ix) &
                               + lambda_w * phas_qp + lambda_x * Iphas_Um 
                 TEB(3,a_ix) = TEB(3,a_ix) &
                              +  lambda_w * phas_Up - lambda_x * Iphas_Qm
                
                 else                

                 TEB(2,a_ix) = TEB(2,a_ix) &
                               + lambda_w * phas_Qm + lambda_x*Iphas_Up
                 TEB(3,a_ix) = TEB(3,a_ix) &
                              +  lambda_w * phas_Um - lambda_x*Iphas_Qp
       
                 end if
              end if ! l allowed by spin or zero

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_rec = h%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_2)  >  ovflow) then
                   lam_0=lam_0/ovflow
                   lam_1=lam_1/ovflow
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = scalefactor(scalem+scalel)*lam_mm/ovflow
                elseif (abs(lam_2)  <  unflow) then
                   lam_0=lam_0*ovflow
                   lam_1=lam_1*ovflow
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = scalefactor(scalem+scalel)*lam_mm/ovflow
                endif

             enddo ! loop on l
          enddo ! loop on m
       endif ! test on cut sky
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixfreerecfac(h)
    deallocate(lam_fact,lfac)
    deallocate(normal_l,twocthlm1)
    deallocate(phas_nq,phas_nu)
    deallocate(phas_sq,phas_su)
    deallocate(phas_n)
    deallocate(phas_s)
    deallocate(ring)
#ifdef MPIPIX
    if (h%MpiId>0) deallocate(map2N,map2S)
    starttime = getetime()
    if (H%MpiId==0) then
     call MPI_REDUCE(MPI_IN_PLACE,TEB,size(TEB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    else
     call MPI_REDUCE(TEB,MPI_IN_PLACE,size(TEB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    end if
    if (debugmsgs>1) print *,code//' time at reduce ', h%MpiId, getetime() -starttime
    if (h%MpiId == 0) call packTEB2TEB(TEB,alm_TEB, nlmax)
    if (debugmsgs>0 .and. h%MpiId==0) print *,code //' time: ',getetime() - initime
#else
    call packTEB2TEB(TEB,alm_TEB, nlmax)
#endif
    deallocate(TEB)

  end subroutine map2polalm


  subroutine polalm2map(H,inlmax, alm_TEB, map_TQU)
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B), INTENT(IN) :: inlmax
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB
    REAL(SP), DIMENSION(:,:), pointer :: map2N,map2S
    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n, b_s

    CHARACTER(LEN=*), PARAMETER :: code = 'POLALM2MAP'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(0:4*H%nside-1) :: ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer a_ix, nalms
#ifdef MPIPIX
    double precision Initime
    integer i
#endif    !=======================================================================

    !     --- allocates space for arrays ---


     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
     call SyncInts(nlmax)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(TEB(3,nalms))
     if (H%MpiId==0) call TEB2PackTEB(alm_TEB,TEB,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
     if (H%MpiId==0) then
      map2N => map_TQU
      map2S => map_TQU
     else
      allocate(map2N(H%North_Start(H%MpiId):H%North_Start(H%MpiId)+H%North_Size(H%MpiId)-1,3),stat = status) 
      if (status /= 0) call die_alloc(code,'map2')   
      allocate(map2S(H%South_Start(H%MpiId):H%South_Start(H%MpiId)+H%South_Size(H%MpiId)-1,3),stat = status) 
      if (status /= 0) call die_alloc(code,'map2')
     end if   
#else
     map2S => map_TQU 
     map2N => map_TQU
#endif


    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   lam_fact = lam_fact * 2 !HealPix polarization def
   

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
       normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo

 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !     --------------------------------------------

    do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       mmax_ring = get_mmax(nlmax,sth) 

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac     !  actual lambda_mm      

          a_ix = a_ix + 1

          b_n = lam_lm * TEB(1,a_ix)
          b_s = b_n

          !l=m special case
          if (m >=2) then
              lambda_w = - 2.0_dp *(normal_l(m) * lam_lm * (fm - fm2) ) * ( one_on_s2 - 0.5_dp )
              lambda_x =  ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              b_n_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_s_Q =  par_lm*(lambda_w * TEB(2,a_ix) - zi_lam_x * TEB(3,a_ix))

              b_n_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
              b_s_U = par_lm*(lambda_w * TEB(3,a_ix) + zi_lam_x * TEB(2,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1

             factor = lam_lm * TEB(1,a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U
          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H,nlmax,b_north,nph,ring,kphi0,mmax_ring)   
       map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,1) = ring(0:nph-1)
       call spinring_synthesis(H,nlmax, b_north_Q, nph, ring, kphi0,mmax_ring)
       map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,2) = ring(0:nph-1)
       call spinring_synthesis(H,nlmax, b_north_U, nph, ring, kphi0,mmax_ring)
       map2N(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,3) = ring(0:nph-1)
  
       if (ith  <  2*nsmax) then
          call spinring_synthesis(H,nlmax, b_south, nph, ring, kphi0,mmax_ring)
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1,1) = ring(0:nph-1)
          call spinring_synthesis(H,nlmax, b_south_Q, nph, ring, kphi0,mmax_ring)
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1,2) = ring(0:nph-1)
          call spinring_synthesis(H,nlmax, b_south_U, nph, ring, kphi0,mmax_ring)
          map2S(H%istart_south(ith):H%istart_south(ith)+nph-1,3) = ring(0:nph-1)
       endif

    enddo    ! loop on cos(theta)


    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)

    deallocate(TEB)

#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiId
    StartTime = Getetime()
    if (H%MpiSize>1) then
     do i=1,3
        if (H%MpiID==0) then
         call MPI_GATHERV(MPI_IN_PLACE,H%North_Size(H%MpiId),SP_MPI, &
          map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
         call MPI_GATHERV(MPI_IN_PLACE,H%South_Size(H%MpiId),SP_MPI, &
          map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
        else
         call MPI_GATHERV(map2N(H%North_Start(H%MpiId),i),H%North_Size(H%MpiId),SP_MPI, &
          map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
         call MPI_GATHERV(map2S(H%South_Start(H%MpiId),i),H%South_Size(H%MpiId),SP_MPI, &
          map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
        end if
     end do
    end if
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code // ' Time: ', GeteTime() - iniTime
    if (H%MpiID>0) deallocate(map2N,map2S) 
#endif

  end subroutine polalm2map


   !=======================================================================
  subroutine scalalm2bispectrum(H, inlmax, inlmax_bi, inL1_max, alm, InSlice)
    !=======================================================================
    ! calculate bispectrum estimator by taking real space products of map rings
    !=======================================================================
    use MPIstuff
    use AMLUtils, only : GetThreeJs   
    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax, inL1_max, inlmax_bi
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
    REAL(DP), INTENT(IN), target :: InSlice(inlmax_bi,inL1_max)
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax, lmax_bi          ! alm related
    INTEGER(I4B) :: nph, kphi0                         ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2
    REAL(DP) :: f2m, corfac
    COMPLEX(DPC) :: factor, factor2
    integer L1_max

    CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2BISPECTRUM'
    COMPLEX(DPC), allocatable :: b_north(:,:), b_flipped(:)
    INTEGER(I4B) :: mmax_ring !,  par_lm
    integer nalms, a_ix
    integer L1,L2,L3, min_l, max_l
    real(DP) aring1(0:4*H%nside-1),aring2(0:4*H%nside-1), tmp
    REAL(DP), pointer :: Slice(:,:)
    REAL(DP), allocatable :: a3j(:)
    REAL(SP), allocatable :: ring(:,:)
    integer, parameter :: L1_min=4, lmin_calc=1200
#ifdef MPIPIX    
    double precision Initime
    integer status
#endif      
!=======================================================================
  
     nsmax = H%nside
     nlmax = inlmax
     L1_Max = inL1_max
     lmax_bi=inlmax_bi

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiId==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call SyncInts(nlmax,L1_Max,lmax_bi)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms))
     if (H%MpiId==0) then
      call Alm2PackAlm(alm,alm2,nlmax)
      Slice=>InSlice
     else
      allocate(Slice(lmax_bi,L1_max))
     end if
     Slice=0
    
#ifdef MPIPIX
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiId, GeteTime() - StartTime
#endif

    call HealpixInitRecfac(H,nlmax)
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    allocate(b_north(0:H%lmax,0:H%lmax))
    allocate(b_flipped(0:H%lmax))
    b_north=0
    allocate(ring(0:4*H%nside-1,0:H%lmax))
    
    do ith =H%ith_start(H%MpiId), H%ith_end(H%MpiId)  
       !        cos(theta) in the pixelisation scheme

       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       mmax_ring = get_mmax(nlmax,sth) 

       lam_mm = sq4pi_inv ! lambda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          f2m = 2.0_dp * m

          !           ---------- l = m ----------
!          par_lm = 1  ! = (-1)^(l+m)
          if (m >= 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)*lam_mm/OVFLOW
  
          lam_lm = corfac
          a_ix = a_ix + 1
          b_north(m,m) = lam_lm * alm2(a_ix)
         
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp 
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec

          do l = m+1, nlmax-1, 2
             
            a_ix = a_ix+1

            lam_0 = lam_1 / a_rec
            lam_1 = lam_2
    
            a_rec = H%recfac(a_ix)
            lam_2 = (cth * lam_1 - lam_0) * a_rec
            
            factor = (lam_1*corfac) * alm2(a_ix)
            factor2 = (lam_2*corfac) * alm2(a_ix+1)

            b_north(l,m) = factor
            b_north(l+1,m) = factor2
            
            lam_0 = lam_1 / a_rec
            lam_1 = lam_2
            a_ix = a_ix+1
            a_rec = H%recfac(a_ix)
            lam_2 = (cth * lam_1 - lam_0) * a_rec
           
             if (abs(lam_1+lam_2) > OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             elseif (abs(lam_1+lam_2) < UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)*lam_mm/OVFLOW
             endif
          enddo
          if (mod(nlmax-m,2)==1) then
                 a_ix = a_ix+1
                 b_north(nlmax,m) = corfac*lam_2*alm2(a_ix)
          end if

       enddo
 
       do L=L1_min,lmax_bi,L1_min
        if (L > L1_max .and. L<lmin_calc) cycle
        b_flipped = b_north(L,:)
        call spinring_synthesis(H,nlmax,b_flipped,nph,ring(0,L),kphi0,min(L,mmax_ring))   ! north hemisph. + equator
       end do
        
       tmp = pi / (3.0_dp * nsmax * nsmax)
       if (ith < 2*nsmax) then
        !Add north and south
          tmp=tmp * 2.d0
       end if         
       do L1=L1_min,L1_max,L1_min
        aring1(0:nph-1) = ring(0:nph-1,L1)         
        do L2=L1,lmax_bi,L1_min
         if (L2<lmin_calc) cycle
         aring2(0:nph-1) = aring1(0:nph-1) *ring(0:nph-1,L2) 
         L3=L2+L1
         if (L3>lmax_bi) cycle
         Slice(L2,L1)=Slice(L2,L1)+ tmp*sum(aring2(0:nph-1)*ring(0:nph-1,L3))
   
!         min_l = max(abs(l1-l2),l2)
!         if (mod(l1+l2+min_l,2)/=0) then
!              min_l = min_l+1
!         end if 
!         max_l = min(nlmax,L1+L2) 
!        
!         do l3=min_l,max_l ,2
!          a_ix=a_ix+1
!          Bispectrum(a_ix,L1)=Bispectrum(a_ix,L1) + sum(aring2(0:nph-1)*ring(0:nph-1,L3))
!         end do !L3
        end do  !L2
       end do  !L1  
 
  
       
    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecFac(H)
    deallocate(ring)
    deallocate(b_north,b_flipped)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code //' Gather ',H%MpiId
    
    StartTime = Getetime()
    if (H%MpiSize>1) then
    
     if (H%MpiId==0) then
      call MPI_REDUCE(MPI_IN_PLACE,Slice,size(Slice),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,l) 
     else
      call MPI_REDUCE(Slice,MPI_IN_PLACE,size(Slice),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,l) 
      deallocate(Slice)
     end if
    
    end if
    
    if (DebugMsgs>1) print *,code //' Done Reduce ',H%MpiId, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiId==0) print *,code //' Time :', GeteTime()-IniTime
#endif
   if (H%MpiID==0) then
      allocate(a3j(0:lmax_bi*2+1))
      do L1=L1_min,L1_max, L1_min
       do L2= L1, lmax_bi, L1_min
         if (L2<lmin_calc) cycle
         L3=L2+L1
         if (L3>lmax_bi) cycle
         call GetThreeJs(a3j(abs(l2-l1)),l1,l2,0,0)
         tmp = real((2*L1+1)*(2*L2+1),dp)*(2*L3+1)/fourpi 
         Slice(L2,L1)=Slice(L2,L1) / a3j(L3)**2/tmp 
        end do
       end do 
      deallocate(a3j) 
   end if
!   
  end subroutine scalalm2bispectrum



 subroutine PackAlm2Alm(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(Out), DIMENSION(1:1,0:nlmax,0:nlmax) :: almout
   COMPLEX(SPC), INTENT(IN), DIMENSION(((nlmax+1)*(nlmax+2))/2) :: almin

   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(1,l,m) = almin(a_ix)
     end do
    end do  

 end  subroutine PackAlm2Alm

 subroutine Alm2PackAlm(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(in), DIMENSION(1:1,0:nlmax,0:nlmax) :: almin
   COMPLEX(SPC), INTENT(out), DIMENSION(((nlmax+1)*(nlmax+2))/2) :: almout
   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(a_ix) = almin(1,l,m) 
     end do
    end do  

 end  subroutine Alm2PackAlm


 subroutine PackEB2EB(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(Out), DIMENSION(1:2,0:nlmax,0:nlmax) :: almout
   COMPLEX(SPC), INTENT(IN), DIMENSION(1:2,((nlmax+1)*(nlmax+2))/2) :: almin

   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(:,l,m) = almin(:,a_ix)
     end do
    end do  

 end  subroutine PackEB2EB

 subroutine EB2PackEB(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(in), DIMENSION(1:2,0:nlmax,0:nlmax) :: almin
   COMPLEX(SPC), INTENT(out), DIMENSION(1:2,((nlmax+1)*(nlmax+2))/2) :: almout
   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(:,a_ix) = almin(:,l,m) 
     end do
    end do  

 end  subroutine EB2PackEB


 subroutine TEB2PackTEB(almin,TEBout,nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(in), DIMENSION(1:3,0:nlmax,0:nlmax) :: almin
   COMPLEX(SPC), INTENT(out), DIMENSION(1:3,((nlmax+1)*(nlmax+2))/2) :: TEBout
  
   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      TEBout(:,a_ix) = almin(:,l,m) 
     end do
    end do  

 end  subroutine TEB2PackTEB


subroutine PackTEB2TEB(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(Out), DIMENSION(1:3,0:nlmax,0:nlmax) :: almout
   COMPLEX(SPC), INTENT(IN), DIMENSION(1:3,((nlmax+1)*(nlmax+2))/2) :: almin

   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(:,l,m) = almin(:,a_ix)
     end do
    end do  

 end  subroutine PackTEB2TEB

  subroutine healpix_sleepMPI(h)
    use mpistuff
    use AMLutils
    type (healpixinfo) :: h
    character(len=*), parameter :: code = 'WAIT'
    
#ifdef MPIPIX
     if (h%MpiId==0) then 
      call sendmessages(h,code)
     else
      call MpiQuietWait   
     end if
#endif
  end  subroutine healpix_sleepMPI

   
   subroutine healpix_wakeMPI
    use AMLutils
    call MpiWakeQuietWait
   end  subroutine healpix_wakeMPI
    
#ifdef MPIPIX

  subroutine MessageLoop(H)
   use MPIStuff
    Type (HealpixInfo) :: H
    character (LEN=64) :: Msg
    REAL(SP),   DIMENSION(1) :: dummymap
    REAL(SP),   DIMENSION(1,3) :: dummymapTQU
    COMPLEX(SPC),   DIMENSION(1) :: dummymapC
    COMPLEX(SPC),   DIMENSION(1,1,1)  :: dummyalm
    REAL(DP), DIMENSION(1,1) :: dummyslice
    Type (HealpixMapArray) :: dummymaps(1)
    Type(HealpixPackedScalAlms) :: dummyscalalms
    Type(HealpixPackedAlms) :: dummyalms     
    integer :: i = 0

     do
      Msg = ''
      call MPI_BCAST(Msg,64,MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)  !***

      if (DebugMsgs>1) print *,'Got message ', H%MpiId, ' '//trim(Msg)

      if (Msg=='EXIT') return
      if (Msg == 'WAIT') then
          call healpix_sleepMPI(H)
      else if (Msg == 'SCALALM2MAP') then
         call scalalm2map(H,i,dummyalm, dummymap)
      else if (Msg == 'ALM2GRADIENTMAP') then
         call alm2GradientMap(H, i , dummyalm, dummymapC)
      else if (Msg == 'SPINALM2MAP') then
        call SpinAlm2Map(H,i,dummyalm, dummymapC, 1)
      else if (Msg == 'MAP2SCALALM') then
        call map2scalalm(H, i, dummymap, dummyalm, -1.d0)   
      else if (Msg=='MAP2SPINALM') then
        call map2spinalm(H, i, dummymapC, dummyalm, 0, -1.d0)   
      else if (Msg=='SCALALM2LENSEDMAP') then
         call scalalm2LensedMap(H,i,dummyalm, dummymapC, dummymap)
      else if (Msg=='ALM2LENSEDMAP') then
         call alm2LensedMap(H,i,dummyalm, dummymapC, dummymapTQU)
      else if (Msg=='ALM2LENSEDMAPINTERP') then
         call alm2LensedmapInterp(H,i,dummyalm, dummymapC, dummymapTQU)
      else if (Msg=='SCALALM2LENSEDMAPINTERP') then
         call scalalm2LensedmapInterp(H,i,dummyalm, dummymapC, dummymap)
      else if (Msg=='SCALALM2LENSEDMAPINTERPCYL') then
         call scalalm2LensedmapInterpCyl(H,i,dummyalm, dummymapC, dummymap)
      else if (Msg=='ALM2LENSEDMAPINTERPCYL') then
         call alm2LensedmapInterpCyl(H,i,dummyalm, dummymapC, dummymapTQU)
      else if (Msg=='MAP2POLALM') then
         call Map2PolAlm(H,i,dummymapTQU, dummyalm, -1.d0)
      else if (Msg=='POLALM2MAP') then
         call polalm2map(H,i,dummyalm, dummymapTQU)
       else if (Msg=='ALM2LENSEDQUADCONTRIB') then
          call alm2LensedQuadContrib(H, i, dummyalm, dummymapC, dummymap)
      else if (Msg=='MAPARRAY2PACKEDSCALALMS') then
          call  maparray2packedscalalms(H,i, dummymaps, dummyscalalms, 1, .false.)
      else if (Msg=='MAPARRAY2PACKEDPOLALMS') then
          call  maparray2packedpolalms(H,i, dummymaps, dummyalms, 1, .false.)
      else if (Msg=='SCALALM2BISPECTRUM') then
          call  scalalm2bispectrum(H, i, i, i, dummyalm, dummyslice)
      end if
     end do 
   
  end subroutine MessageLoop
 
  subroutine SendMessages(H, MsgIn)
   use MPIStuff
   Type (HealpixInfo) :: H
   CHARACTER(LEN=*), intent(in) :: MsgIn
   CHARACTER(LEN=64) :: Msg
    
    Msg = MsgIn
    if (DebugMsgs>1) print *,'Send messages '//trim(Msg)
    call MPI_BCAST(Msg,64,MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr) !***
   ! same time as BCAST in MessaegLoop

  end  subroutine SendMessages

#endif
end module spinalm_tools
