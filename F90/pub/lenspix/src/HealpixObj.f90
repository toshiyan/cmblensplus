module HealpixObj
 use healpix_types, ONLY: SP,DP,I4B,SPC
 USE head_fits, ONLY : add_card, get_card
 USE pix_tools, ONLY :  npix2nside, nside2npix, query_disc
 USE fitstools, ONLY : getsize_fits, input_map, read_par, read_dbintab, write_asctab, &
   dump_alms,write_bintab
 USE spinalm_tools
 use AMLutils
 implicit none

 real(SP), parameter :: mK = 1. !internal units/microK
 
 REAL(SP) :: fmissval = -1.6375e-30
 real(DP), parameter :: HO_PI = 3.14159265358979323846264338328d0, &
      HO_twopi=2*HO_pi, HO_fourpi=4*HO_pi

 integer, parameter :: ord_ring = 1, ord_nest = 2

 integer, parameter ::  C_T = 1, C_E = 2, C_B = 3, C_C = 4
 integer, parameter :: C_TB = 5, C_EB = 6 !namikawa

 logical :: RandInited = .false.

 integer, parameter :: nospinmap = 0

 type HealpixMap

  integer(I4B) npix, nmaps, ordering, nside, type
  REAL(SP), DIMENSION(:,:),   POINTER :: TQU  => NULL() 
  COMPLEX(SP), dimension(:), pointer :: SpinField => NULL() 
!  REAL(SP), DIMENSION(:), pointer :: SpinQ, SpinU
  REAL(SP), dimension(:), pointer :: Phi => NULL() 
  integer :: spin
  logical :: HasPhi
 end type HealpixMap


 type HealpixAlm

   integer(I4B) lmax, npol, spin
   logical HasPhi
   COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: TEB  => NULL()    !T,E,B index, l,m 
   COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: SpinEB  => NULL()  ! E,B index, l, m
   COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: Phi => NULL()  !Has to be 3D array for Healpix defs

 end type HealpixAlm

 type HealpixPower
   !Raw Cls
   !read from text files in l(l+1)C_l/(2pi microK^2)
   !PhiCl is lensing potential phi-phi and phi-T and phi-E. 
   !Phi is read in above units, but stored here as dimensionless
   integer(I4B) lmax
   logical pol, lens
   REAL(SP), DIMENSION(:,:),   POINTER :: Cl => NULL() 
   REAL(SP), DIMENSION(:,:),  POINTER :: PhiCl   => NULL() 

 end  type HealpixPower
 
 type HaarComponents
   
   integer order
   logical pol
   type(HealpixMap) degraded
   type(HealpixMap), dimension(:), pointer :: details  => NULL() 
  
 end  type HaarComponents

contains

   subroutine HealpixPower_Init(P, almax, pol, dolens, nofree)
    Type(HealpixPower) P
    integer, intent(in) :: almax
    logical, intent(in) :: pol
    logical, intent(in), optional :: dolens
    logical, intent(in), optional :: nofree
    logical dofree
    
    if (present(nofree)) then
     dofree=.not. nofree
    else
     dofree = .true. 
    end if  
      
     if (dofree) call HealpixPower_Free(P)
     P%lmax = almax
     P%pol = pol
     P%lens = .false.
     if (present(dolens)) then
       P%lens = dolens
     end if
     if (pol) then
       !allocate(P%Cl(0:almax,4))  
       allocate(P%Cl(0:almax,6)) !namikawa
     else
       allocate(P%Cl(0:almax,1))  
     end if
     P%Cl=0
     if (P%lens) then
      allocate(P%PhiCl(0:almax,3))      
      P%PhiCl = 0
     end if

   end subroutine HealpixPower_Init

   subroutine HealpixPower_Nullify(P) 
    Type(HealpixPower) P

    nullify(P%Cl, P%PhiCl)
   
   end subroutine HealpixPower_Nullify
   
   subroutine HealpixPower_Free(P)
    Type(HealpixPower) P
    integer status

    deallocate(P%Cl, stat = status)
    deallocate(P%PhiCl, stat = status)
    call HealpixPower_Nullify(P)

   end subroutine HealpixPower_Free

   subroutine HealpixPower_Assign(P, Pin)
    Type(HealpixPower) P, Pin

    call HealpixPower_Init(P, Pin%lmax, Pin%pol, Pin%lens)
    P%Cl = Pin%Cl
    if (Pin%lens) P%PhiCl = Pin%PhiCl  

   end subroutine HealpixPower_Assign


subroutine HealpixPower_ReadFromTextFile(P, f, lmax, pol, dolens)
  use AMLutils
  Type(HealpixPower) P
  character(LEN=*), intent(IN) :: f
  integer, intent(in) :: lmax
  logical,intent(in), optional :: dolens, pol
  integer i,l, ColNum
  real(DP) T,E,B,TE, scal, phi, phiT, phiE
  logical tensors, dolens2, pol2
  integer, parameter :: io_unit= 1
  character(LEN=200) :: InLine
  real(sp) test(9)
  real(sp), parameter :: COBE_CMBTemp = 2.726

  open(unit=io_unit,file=f,form='formatted',status='old')

  call HealpixPower_Free(P)
  dolens2=.false.
  if (present(dolens)) then
    dolens2 = dolens
  end if
  pol2 = .false.
  if (present(pol)) then
    pol2 = pol
  end if
  call HealpixPower_Init(P, lmax, pol = pol2, dolens = dolens2)
 
  !See how many columns - includes magnetic polarization if four columns
  read(io_unit,'(a)') InLine
  do l=9,2,-1
    Colnum=l 
    read(InLine,*, end=110) test(1:l)
    exit
110 cycle
  end do 
  tensors = colnum==5 .or. colnum > 7
  if (P%lens .and. (colnum < 6 .or. colnum>8)) call MpiStop('can''t indentify text phi cl columns')
  if (pol2 .and. colnum<3) call mpiStop('No polarization in C_l file' ) 
  rewind io_unit
  B=0
  phiE=0
  P%Cl = 0
  if (dolens2) P%PhiCl = 0
  do i=0,lmax
    if (Colnum==2) then
      read (io_unit,*,end=118) l, T
    else
      if (tensors) then
        if (P%lens) then
          if (colnum ==8) then
            read (io_unit,*,end=118) l, T, E, B , TE, phi, phiT, phiE
          else
            read (io_unit,*,end=118) l, T, E, B , TE, phi, phiT
          end if
        else
          read (io_unit,*,end=118) l, T, E, B , TE
        end if
      else
        if (P%lens) then
          read (io_unit,*,end=118) l, T, E, TE, phi, phiT, phiE
!          read (io_unit,*,end=118) l, T, TE, E, phi, phiT !namikawa
        else
          read (io_unit,*,end=118) l, T, E, TE
!          read (io_unit,*,end=118) l, T, TE, E !namikawa
        end if
      end if
    end if
    if (l==0) then
      P%Cl(l,1) = T
    else if (l<= lmax) then
      scal = twopi/(l*(l+1))
      P%Cl(l,1) = T*scal
      if (pol2 .and. l>=2 ) then
        P%Cl(l,2) = E*scal
        P%Cl(l,3) = B*scal
        P%Cl(l,4) = TE*scal
      end if
      if (P%Lens .and. l>=1) then
        P%PhiCl(l,1) = phi * twopi/real(l*(l+1),dp) !namikawa
        P%PhiCl(l,2) = phiT * twopi/real(l*(l+1),dp) !namikawa
        P%PhiCl(l,3) = phiE * twopi/real(l*(l+1),dp) !namikawa
      end if
    end if
  end do
  118 close(io_unit)

  !open(unit=20,file="confirm-input.dat",status="replace")
  !do l = 2, lmax
  !  write(20,"(I6,2X,6(E14.7,2X))") l, P%Cl(l,1), P%Cl(l,2), P%Cl(l,3), P%Cl(l,4), P%PhiCl(l,1), P%PhiCl(l,2)
  !end do
  !close(20)


end subroutine HealpixPower_ReadFromTextFile

subroutine HealpixPower_Write(P,fname) 
  use AMLutils
  Type(HealpixPower) P
  character(Len=*), intent(in) :: fname
  integer l
  real fac

  call CreateTxtFile(fname,1)
  if (P%pol) then
    if (P%lens) then
      write (1,'(1I7,7E17.7)') 0,P%Cl(0,:) *mK**2, 0.,0.,0. 
    else
      write (1,'(1I7,4E17.7)') 0,P%Cl(0,:) *mK**2
    end if 
  else 
    if (P%lens) then
      write (1,'(1I7,3E17.7)') 0,P%Cl(0,1) *mK**2, 0.,0.
    else
      write (1,'(1I7,1E17.7)') 0,P%Cl(0,1) *mK**2         
    end if
  end if

  do l=1,P%lmax
    fac = l*(l+1)
    if (P%pol) then
      if (P%lens) then
!        write (1,'(1I7,7E17.7)') l,fac*P%Cl(l,:)/twopi *mK**2, &
!              fac**2*P%PhiCl(l,1)/twopi,fac**1.5*P%PhiCl(l,2:3)/twopi * mK 
        write (1,'(1I7,7E17.7)') l,fac*P%Cl(l,:)/twopi *mK**2, &
          fac*P%PhiCl(l,1)/twopi,fac*P%PhiCl(l,2:3)/twopi * mK !namikawa
      else
        write (1,'(1I7,4E17.7)') l,fac*P%Cl(l,:)/twopi *mK**2
      end if
    else
      if (P%lens) then
!        write (1,'(1I7,3E17.7)') l,fac*P%Cl(l,1)/twopi *mK**2, &
!          fac**2*P%PhiCl(l,1)/twopi,fac**1.5*P%PhiCl(l,2)/twopi * mK 
        write (1,'(1I7,3E17.7)') l,fac*P%Cl(l,1)/twopi *mK**2, &
          fac*P%PhiCl(l,1)/twopi,fac*P%PhiCl(l,2)/twopi * mK !namikawa
      else
        write (1,'(1I7,1E17.7)') l,fac*P%Cl(l,1)/twopi *mK**2
      end if
    end if
  end do
  close(1)

end subroutine HealpixPower_Write

subroutine HealpixPower_AddPower(Ptotal, P, AddPhi)
  !Adds P to PTotal (useful for getting means over realisations). 
  Type(HealpixPower) P, Ptotal
  logical, intent(in) :: AddPhi
      
  if (all(shape(P%Cl) /= shape(Ptotal%Cl))) &
    call MpiStop('HealpixPower_AddPower: must have same sized power spectra') 
    Ptotal%Cl = Ptotal%Cl + P%Cl
  if (AddPhi) then
    if (.not. PTotal%lens  .or. .not. P%lens) call MpiStop('HealpixPower_AddPower: must both have phi')
    Ptotal%PhiCl = PTotal%phiCl + P%PhiCl 
  end if

end subroutine HealpixPower_AddPower

subroutine HealpixPower_Smooth(P,fwhm, sgn)
  Type(Healpixpower) :: P
  integer l, sn
  integer, intent(in), optional :: sgn
  real(dp) xlc,sigma2,fwhm

  if (present(sgn)) then
    sn = sgn
  else
    sn = -1
  end if
   
  xlc= 180*sqrt(8.*log(2.))/HO_pi
  sigma2 = (fwhm/xlc)**2

  do l=2,P%lmax
    P%Cl(l,:) =  P%Cl(l,:)*exp(sn*l*(l+1)*sigma2)
  end do

end subroutine HealpixPower_Smooth

subroutine HealpixPower_Smooth_Beam(P,beam, sgn)
  Type(Healpixpower) :: P
  integer l, sn
  integer, intent(in), optional :: sgn
  real(dp), intent(in) :: beam(0:)

  if (present(sgn)) then
    sn = sgn
  else
    sn = -1
  end if
   
  do l=0,P%lmax
    if (sn==-1) then
      P%Cl(l,:) =  P%Cl(l,:)* beam(l)**2
    else
      P%Cl(l,:) =  P%Cl(l,:)/ beam(l)**2
    end if
  end do

end subroutine HealpixPower_Smooth_Beam

subroutine HealpixPower_Smooth_Beam2(P,beam1,beam2, sgn)
  Type(Healpixpower) :: P
  integer l, sn
  integer, intent(in), optional :: sgn
  real(dp), intent(in) :: beam1(0:), beam2(0:)

  if (present(sgn)) then
    sn = sgn
  else
    sn = -1
  end if
   
  do l=0,P%lmax
    if (sn==-1) then
      P%Cl(l,:) =  P%Cl(l,:)* beam1(l)*beam2(l)
    else
      P%Cl(l,:) =  P%Cl(l,:)/ ( beam1(l)*beam2(l))
    end if
  end do

end subroutine HealpixPower_Smooth_Beam2

subroutine HealpixAlm2Power(A,P)
  Type(HealpixPower) P
  Type(HealpixAlm) :: A
  integer l,i,ix

  call HealpixPower_Init(P,A%lmax,A%npol==3, A%HasPhi)
    
  P%Cl(0:1,:) = 0
  do l=0, P%lmax
    if (l<2) then
      ix= 1
    else
      ix = A%npol
    end if
    do i = 1, ix
      P%Cl(l,i) = ( REAL(A%TEB(i,l,0))**2 &
        + 2.*SUM(A%TEB(i,l,1:l)*CONJG(A%TEB(i,l,1:l)) )) / (2.*l + 1.)
!      write(*,*) l, A%TEB(2,l,l), CONJG(A%TEB(i,l,l)) !namikawa
    end do
    if (ix==3) then
      P%Cl(l,4) = ( REAL(A%TEB(1,l,0))*REAL(A%TEB(2,l,0)) &
        + 2.*SUM(real(A%TEB(1,l,1:l)*CONJG(A%TEB(2,l,1:l))) ) &
          ) / (2.*l + 1.)
    end if
  end do 
  if (A%HasPhi) then
    do l=0, P%lmax
      P%PhiCl(l,1) = ( REAL(A%Phi(1,l,0))**2 &
        + 2.*SUM(A%Phi(1,l,1:l)*CONJG(A%Phi(1,l,1:l)) )) / (2*l + 1)
      !T-phi
      P%PhiCl(l,2) = ( REAL(A%TEB(1,l,0))*REAL(A%Phi(1,l,0)) &
        + 2.*SUM(real(A%TEB(1,l,1:l)*CONJG(A%Phi(1,l,1:l))) )) /(2*l + 1)
      if (l>1 .and. A%npol >2) then
        !E-phi
        P%PhiCl(l,3) = ( REAL(A%TEB(2,l,0))*REAL(A%Phi(1,l,0)) &
          + 2.*SUM(real(A%TEB(2,l,1:l)*CONJG(A%Phi(1,l,1:l))) )) /(2*l + 1)
      end if 
    end do 
  end if
end subroutine HealpixAlm2Power
   
subroutine HealpixAlm2CrossPower(A,A2, P)
    Type(HealpixPower) P
    Type(HealpixAlm) :: A, A2
    integer l,i,ix

    if (A%lmax /= A2%lmax) call MpiStop('HealpixAlm2CrossPower: mismatched lmax')
    if (A%npol /= A2%npol) call MpiStop('HealpixAlm2CrossPower: different pol content')
    call HealpixPower_Init(P,A%lmax,A%npol==3)
    
    P%Cl(0:1,:) = 0
    do l=0, P%lmax
     if (l<2) then
      ix= 1
     else
      ix = A%npol
     end if
     do i = 1, ix
      P%Cl(l,i) = ( REAL(A%TEB(i,l,0)*A2%TEB(i,l,0)) &
            + 2.*SUM(A%TEB(i,l,1:l)*CONJG(A2%TEB(i,l,1:l)) )) / (2.*l + 1.)
     end do
     if (ix==3) then
        P%Cl(l,4) = ( REAL(A%TEB(1,l,0))*REAL(A2%TEB(2,l,0)) &
            + 2.*SUM(real(A%TEB(1,l,1:l)*CONJG(A2%TEB(2,l,1:l))) ) &
           ) / (2.*l + 1.)
      end if
    end do 

   end subroutine HealpixAlm2CrossPower


   subroutine HealpixAlm_Init(A,almax,npol,spinmap, HasPhi)
     Type(HealpixAlm) :: A
     integer, intent(in) :: almax
     integer, intent(in), optional :: npol, spinmap
     logical, intent(in), optional :: HasPhi
     integer status

     call HealpixAlm_Free(A)

     A%lmax = almax
     if (present(npol)) then
      A%npol = npol
     else
      A%npol = 1
     end if

     if (A%npol /= 0) then 
      ALLOCATE(A%TEB(1:A%npol, 0:almax, 0:almax),stat = status)
      if (status /= 0) call MpiStop('No Mem: HealpixAlm_Init lmax = '//IntToStr(almax))
      A%TEB=0
     end if

     if (present(spinmap)) then
         if (spinmap /= nospinmap) then
          if (spinmap<1 .or. spinmap > 3) call mpiStop( 'Spin must be 0<spin<4')
          ALLOCATE(A%SpinEB(2, 0:almax, 0:almax),stat = status)
          if (status /= 0) call MpiStop('No Mem: HealpixAlm_Init spinmap')
         end if
         A%spin= spinmap     
     else
      A%spin = nospinmap      
     end if

     if (present(HasPhi)) then
       A%HasPhi = HasPhi
     else
       A%HasPhi = .false.   
     end if

     if (A%HasPhi) then
          ALLOCATE(A%Phi(1:1,0:almax, 0:almax),stat = status)
          A%Phi = 0
          if (status /= 0) call MpiStop('No Mem: HealpixAlm_Init phi')
     end if

   end subroutine HealpixAlm_Init


  subroutine HealpixAlm_Assign(AOut, Ain, max_pol)
   Type(HealpixAlm) :: AOut, Ain
   integer, intent(in), optional :: max_pol
   integer status
   integer maxp

   if (present(max_pol)) then
    maxp=max_pol
   else
    maxp =Ain%npol
   end if 
   call HealpixAlm_Free(AOut)
   Aout = Ain 
   nullify(AOut%TEB, AOut%SpinEB, AOut%Phi)
   if (maxp>0) then
      ALLOCATE(AOut%TEB(1:maxp, 0:AOut%lmax, 0:AOut%lmax),stat = status)
      if (status /= 0) call MpiStop('No Mem: HealpixAlm_Assign')
     AOut%TEB(1:maxp,:,:)= Ain%TEB(1:maxp,:,:)
   end if
   if (AIn%spin /= nospinmap) then
        ALLOCATE(AOut%SpinEB(2, 0:AOut%lmax, 0:AOut%lmax),stat = status)
        if (status /= 0) call MpiStop('No Mem: HealpixAlm_Assign')
        AOut%SpinEB = Ain%SpinEB
   end if
   if (AIn%HasPhi) then
       ALLOCATE(AOut%Phi(1:1,0:AOut%lmax, 0:AOut%lmax),stat = status)
      if (status /= 0) call MpiStop('No Mem: HealpixAlm_Assign')
      AOut%Phi = Ain%Phi   
   end if 
  
  end subroutine HealpixAlm_Assign

  subroutine HealpixAlm_Nullify(A)
   Type(HealpixAlm) :: A

     nullify(A%TEB)
     nullify(A%SpinEB)
     nullify(A%Phi)
 
  end subroutine HealpixAlm_Nullify


  subroutine HealpixAlm_Free(A)
   Type(HealpixAlm) :: A
   integer status

     deallocate(A%TEB,stat=status)
     deallocate(A%SpinEB,stat=status)
     deallocate(A%Phi,stat=status)
     call HealpixAlm_Nullify(A)
     
  end subroutine HealpixAlm_Free

  subroutine HealpixAlm_PhiOnly(A)
   Type(HealpixAlm) :: A
   integer status

     deallocate(A%TEB,stat=status)
     deallocate(A%SpinEB,stat=status)
     nullify(A%TEB)
     nullify(A%SpinEB)
     A%spin = nospinmap      
     A%npol = 0 
  end subroutine 


  subroutine HealpixAlm_GradientOf(A, B, field, updown)
   type(HealpixAlm) :: A, B
   integer l
   character(LEN=*), intent(in) :: field
   character(LEN=*), intent(in), optional :: updown
   logical Div  
   integer spin
   
   if (field(1:1) /= 'S') then
    spin = 1
   else
    if (.not. present(updown))  call MpiStop('HealpixAlm_GradientOf: Must say which derivative')
    Div = updown(1:1) == 'D'

    if (Div) then
  !    spin = A%spin-1
      spin = 1
    else
!     spin = A%spin+1
      spin = 3
    end if
   end if


   call HealpixAlm_Init(B,A%lmax,0,spinmap= spin)

   B%SpinEB = 0
    do l=B%spin, A%lmax
      if (field(1:1)=='P') then
       B%SpinEB(1,l,0:l) = - EB_sign*sqrt(real(l*(l+1),dp))*A%Phi(1,l,0:l)
      else if (field(1:1)=='T') then
       B%SpinEB(1,l,0:l) = - EB_sign*sqrt(real(l*(l+1),dp))*A%TEB(1,l,0:l)
      else if (field(1:1)=='S') then
         if (Div) then
        !Divergence
          B%SpinEB(:,l,0:l) =  sqrt(real((l+2)*(l-1),dp))*A%TEB(2:3,l,0:l)
          else
          !STF of outer product
         B%SpinEB(:,l,0:l) =  -sqrt(real((l+3)*(l-2),dp))*A%TEB(2:3,l,0:l)
         end if
      else
       call MpiStop('HealpixAlm_GradientOf: Unknown field')
      end if 
    end do

  end subroutine HealpixAlm_GradientOf

  subroutine HealpixAlm_PolToSpin2(A)
   Type(HealpixAlm) :: A
   integer status

   deallocate(A%SpinEB, stat=status)
   ALLOCATE(A%SpinEB(2, 0:A%lmax, 0:A%lmax),stat = status)
   A%SpinEB(1:2,:,:) = A%TEB(2:3,:,:)
   A%spin = 2
  end   subroutine HealpixAlm_PolToSpin2

  subroutine HealpixAlm_Spin2ToPol(A)
   Type(HealpixAlm) :: A
   integer :: status

   if (A%npol==0) then
      A%npol = 3
      ALLOCATE(A%TEB(1:A%npol, 0:A%lmax, 0:A%lmax),stat = status)
      if (status /= 0) call MpiStop('No Mem: HealpixAlm_Spin2ToPol')
      A%TEB=0
   end if
   A%TEB(2:3,:,:) = A%SpinEB(1:2,:,:)

  end subroutine HealpixAlm_Spin2ToPol



  subroutine HealpixAlm_Smooth(A,fwhm, sgn)
   Type(HealpixAlm) :: A
   integer l, sn
   integer, intent(in), optional :: sgn
   real(dp) xlc,sigma2,fwhm
   
   if (present(sgn)) then
     sn = sgn
   else
     sn = -1
   end if
   xlc= 180*sqrt(8.*log(2.))/HO_pi
   sigma2 = (fwhm/xlc)**2

   do l=2,A%lmax
     A%TEB(:,l,:) =  A%TEB(:,l,:)*exp(sn*l*(l+1)*sigma2/2)
   end do

  end subroutine HealpixAlm_Smooth

  subroutine HealpixAlm_Smooth_Beam(A,Beam, sgn)
   Type(HealpixAlm) :: A
   integer l, sn
   real(dp) :: beam(0:)
   integer, intent(in), optional :: sgn
   
   if (present(sgn)) then
     sn = sgn
   else
     sn = -1
   end if

   do l=0,A%lmax
     if (sn ==-1) then
      A%TEB(:,l,:) =  A%TEB(:,l,:)*beam(l)
     else
      A%TEB(:,l,:) =  A%TEB(:,l,:)/beam(l)
     end if
   end do

  end subroutine HealpixAlm_Smooth_Beam


  subroutine HealpixMap_GetAzimCut(M, npix,rad, theta,phi)
   !1 inside disc radius rad centred at theta, phi (radians)
    use pix_tools
    Type(HealpixMap) :: M
    real(dp), intent(in) :: rad, theta, phi
    integer, intent(in) :: npix
    real(dp) vec(3)
    
    call ang2vec(theta,phi,vec)
    call HealpixMap_GetAzimCutVec(M, npix, rad, vec)
    
  end subroutine HealpixMap_GetAzimCut

  subroutine HealpixMap_GetAzimCutVec(M, npix,rad, Vec)
   !inside disc radius rad centred at vec
    Type(HealpixMap) :: M
    real(dp), intent(in) :: rad
    integer, intent(in) :: npix
    integer, dimension(:), allocatable :: listpix
    real(dp), intent(in) :: vec(3)
    integer nlist
  
    call HealpixMap_Init(M,npix,1)
    allocate(listpix(0:npix-1))
    call query_disc(M%nside,vec, rad, listpix,nlist)
    M%TQU = 0
    M%TQU(listpix(0:nlist-1),1) = 1
    deallocate(listpix)     

  end subroutine HealpixMap_GetAzimCutVec


  function HealpixMap_Vec2pix(M, vec) result(pix)
    use pix_tools
    Type(HealpixMap) :: M 
    integer :: pix
    real(dp), intent(in) :: vec(3)
    
    if (M%Ordering == ord_ring) then
     call vec2pix_ring(M%nside,vec,pix) 
    else
     call vec2pix_nest(M%nside,vec,pix) 
    end if

  end function HealpixMap_Vec2pix


  subroutine HealpixMap_Pix2Vec(M, pix,vec)
    use pix_tools
    Type(HealpixMap) :: M
 
    integer, intent(in) :: pix
    real(dp), intent(out) :: vec(3)
    
    if (M%Ordering == ord_ring) then
     call pix2vec_ring(M%nside,pix, vec) 
    else
     call pix2vec_nest(M%nside,pix, vec) 
    end if

  end subroutine HealpixMap_Pix2Vec

  subroutine HealpixMap_Pix2Vertex(M, pix,vertex)
    use pix_tools
    Type(HealpixMap) :: M
 
    integer, intent(in) :: pix
    real(dp), intent(out) :: vertex(3,4)
    real(dp) :: vec(3)
    
    if (M%Ordering == ord_ring) then
     call pix2vec_ring(M%nside,pix, vec, vertex) 
    else
     call pix2vec_nest(M%nside,pix, vec, vertex) 
    end if

  end subroutine HealpixMap_Pix2Vertex


  function HealpixMap_Ang2Pix(M, theta, phi) result(pix)
    use pix_tools
    Type(HealpixMap) :: M
    real(dp), intent(in):: theta, phi
    integer pix
    
    if (M%Ordering == ord_ring) then
     call Ang2Pix_ring(M%nside,theta,phi,pix) 
    else
     call Ang2Pix_nest(M%nside,theta,phi,pix) 
    end if

  end function HealpixMap_Ang2Pix
  
  subroutine HealpixMap_Pix2Ang(M, pix, theta, phi) 
    use pix_tools
    Type(HealpixMap) :: M
    integer, intent(in) :: pix
    real(dp), intent(out):: theta, phi
    
    if (M%Ordering == ord_ring) then
     call Pix2Ang_ring(M%nside,pix,theta,phi) 
    else
     call Pix2Ang_nest(M%nside,pix,theta,phi) 
    end if

  end subroutine HealpixMap_Pix2Ang
  

  subroutine Healpix_GetRotation(R, theta, phi, chi)
!Rotate phi about z axis, then rotation by theta about new y axis, then chi about new z axis
    real(dp), intent(in) :: theta, phi, chi 
    real(dp), intent(inout) :: R(3,3)
  
    call MpiStop('Healpix_GetRotation: You''ll have to check this routine')

    R(1,1) = cos(phi)*cos(theta)*cos(chi) - sin(phi)*sin(chi)
    R(1,2) = sin(phi)*cos(theta)*cos(chi) + cos(phi)*sin(chi)
    R(1,3) = -sin(theta)*cos(chi)
    R(2,1) = -cos(phi)*cos(theta)*sin(chi) - sin(phi)*cos(chi)
    R(2,2) = -sin(phi)*cos(theta)*sin(chi) + cos(phi)*cos(chi)
    R(2,3) = sin(phi)*sin(chi)
    R(3,1) = cos(phi)*sin(theta)
    R(3,2) = sin(phi)*sin(theta)
    R(3,3) = cos(theta)

  end subroutine Healpix_GetRotation

 
  subroutine HealpixMap_Rotate(M, MR, theta, phi, chi)
!This is very crude for temp
   use pix_tools
    Type(HealpixMap) :: M, MR
    real(dp), intent(in) :: theta, phi, chi 
    integer i, ix
    real(dp) vec(3), R(3,3)

    call MpiStop('Don''t use this')
    call Healpix_GetRotation(R, theta, phi, chi)
    call HealpixMap_Init(MR, M%npix, M%nmaps)
    call HealpixMap_ForceRing(M)
    do i=0, MR%npix -1
      call pix2vec_ring(MR%nside, i, vec)
      vec = matmul(R,vec)
      call vec2pix_ring(MR%nside, vec, ix)
      MR%TQU(i,:) = M%TQU(ix,:)
    end do

  end subroutine HealpixMap_Rotate


  function HealpixMap_EclipticPixel(M, theta,phi) result (res)
   !Pixel on galactic map for point at given ecliptic coordinates
    use pix_tools
    use coord_v_convert
    Type(HealpixMap) :: M
    real(dp), intent(in) :: theta, phi
    integer res
    real(dp) vec(3), vecout(3)
    
    call ang2vec(theta,phi,vec)
    call xcc_DP_E_TO_G(vec,2000.d0,vecout)
    res = HealpixMap_Vec2pix(M, vecout) 
    
  end function HealpixMap_EclipticPixel


  subroutine HealpixMap_MarkEclipticPlane(M, val)
   !Mark pixels on ecliptic plan by value val, assuming M galactic
    Type(HealpixMap) :: M
    real(dp) :: val, theta, phi, delta
    integer i
    
    theta = HO_pi/2
    delta = 2*HO_pi/ (M%nside*8) 
    do i=1, M%nside*8
     phi = i*delta
     M%TQU(HealpixMap_EclipticPixel(M,theta,phi),:) = val
    end do
    
  end subroutine HealpixMap_MarkEclipticPlane


  subroutine HealpixMap_GalacticToEcliptic(M)
   !Not at all optimal
    Type(HealpixMap) M, M2
    integer i
    real(dp) theta,phi
    
    call HealpixMap_Assign(M2,M)
    do i=0, M%npix-1
      call HealpixMap_Pix2Ang(M, i, theta, phi) 
      M%TQU(i,:) = M2%TQU(HealpixMap_EclipticPixel(M,theta,phi),:)
    end do
    call HealpixMap_Free(M2)
    
  end subroutine HealpixMap_GalacticToEcliptic
 
  subroutine HealpixMap_AddWhiteNoise(M, N_T, N_QU )
  !N_T and N_QU are the N_l of the noise
    use Random
    Type(HealpixMap) :: M
    real(sp), intent(in) :: N_T
    real(sp), intent(in), optional :: N_QU
    real(sp) amp
    integer i
         
    amp = sqrt(N_T*M%npix/(HO_fourpi))
    do i=0, M%npix-1
     M%TQU(i,1)= M%TQU(i,1) + Gaussian1()*amp
    end do
  
    if (present(N_QU) .and. M%nmaps>1) then
        if (M%nmaps /= 3) call MpiStop('HealpixMap_AddWhiteNoise: No polarization in map')
        amp = sqrt(N_QU*M%npix/HO_fourpi)
        do i=0, M%npix-1
         M%TQU(i,2)= M%TQU(i,2) + Gaussian1()*amp
        end do
        do i=0, M%npix-1
         M%TQU(i,3)= M%TQU(i,3) + Gaussian1()*amp
        end do
    end if
    
  end  subroutine HealpixMap_AddWhiteNoise

  subroutine HealpixMap_AddUncorrelatedNoise(M, NoiseMap)
    use Random
    Type(HealpixMap) :: M, NoiseMap
    integer i
    real var
         
    do i=0, M%npix-1
     var= NoiseMap%TQU(i,1)
     if (var>0) M%TQU(i,1)= M%TQU(i,1) + Gaussian1()*sqrt(var)
    end do
  
    if (M%nmaps>1) then
        if (M%nmaps /= 3) call MpiStop('HealpixMap_AddUncorrelatedNoise: No polarization in map')
        if (NoiseMap%nmaps /= 3) call MpiStop('HealpixMap_AddUncorrelatedNoise: No polarization in noise map')
        do i=0, M%npix-1
         var= NoiseMap%TQU(i,2)
         if (var > 0) M%TQU(i,2)= M%TQU(i,2) + Gaussian1()*sqrt(var)
        end do
        do i=0, M%npix-1
         var= NoiseMap%TQU(i,3)
         if (var>0) M%TQU(i,3)= M%TQU(i,3) + Gaussian1()*sqrt(var)
        end do
    end if
    
  end  subroutine HealpixMap_AddUncorrelatedNoise


subroutine HealpixAlm_Sim_TEB(A, P, seed) !namikawa
  use random
  use alm_tools
  use ran_tools
  Type(HealpixAlm) :: A
  Type(HealpixPower) :: P
  integer, intent(in), optional :: seed
  logical :: wantphi
  integer :: wantpol, l, m
  real(dp) :: sqrt2, A11, A21, A22, A31, A32, A33, A41, A42, A43, A44
  complex(dp) :: g1, g2, g3, g4

  if (present(seed)) then
    call InitRandom(seed)
  else
    if (.not. RandInited) call InitRandom
    RandInited = .true.
  end if

  wantpol = 3
  wantphi = .true.

  call HealpixAlm_Init(A,P%lmax, wantpol, HasPhi=wantphi)
  sqrt2 = sqrt(2.)
  A%TEB = 0d0
  A%Phi = 0d0

  do l = 2, P%lmax
    A11 = 0d0;  A21=0d0; A22=0d0; A31=0d0; A32=0d0; A33=0d0; A41=0d0; A42=0d0; A43=0d0; A44=0d0
    g1  = Gaussian1()
    g2  = Gaussian1()
    g3  = Gaussian1()
    g4  = Gaussian1()
    A11 = sqrt(P%cl(l,C_T))
    if (A11/=0) A21 = P%cl(l,C_C)/A11
    A22 = sqrt(P%cl(l,C_E)-A21**2)
    if (A11/=0) A31 = P%cl(l,C_TB)/A11
    if (A22/=0) A32 = (P%cl(l,C_EB)-A21*A31)/A22
    A33 = sqrt(P%cl(l,C_B)-A31**2-A32**2)
    if (A11/=0) A41 = P%phicl(l,2)/A11
    if (A22/=0) A42 = (P%phicl(l,3)-A21*A41)/A22
    A43 = 0d0
    A44 = sqrt(P%phicl(l,1)-A41**2-A42**2-A43**2)
    A%TEB(1,l,0) = g1*A11
    A%TEB(2,l,0) = g1*A21 + g2*A22
    A%TEB(3,l,0) = g1*A31 + g2*A32 + g3*A33
    A%Phi(1,l,0) = g1*A41 + g2*A42 + g3*A43 + g4*A44
    A11 = 0d0;  A21=0d0; A22=0d0; A31=0d0; A32=0d0; A33=0d0; A41=0d0; A42=0d0; A43=0d0; A44=0d0
    do m = 1, l
      g1  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g2  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g3  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      g4  = cmplx(Gaussian1(),Gaussian1())/sqrt2
      A11 = sqrt(P%cl(l,C_T))
      if (A11/=0) A21 = P%cl(l,C_C)/A11
      A22 = sqrt(P%cl(l,C_E)-A21**2)
      if (A11/=0) A31 = P%cl(l,C_TB)/A11
      if (A22/=0) A32 = (P%cl(l,C_EB)-A21*A31)/A22
      A33 = sqrt(P%cl(l,C_B)-A31**2-A32**2)
      if (A11/=0) A41 = P%phicl(l,2)/A11
      if (A22/=0) A42 = (P%phicl(l,3)-A21*A41)/A22
      A43 = 0d0
      A44 = sqrt(P%phicl(l,1)-A41**2-A42**2-A43**2)
      A%TEB(1,l,m) = g1*A11
      A%TEB(2,l,m) = g1*A21 + g2*A22
      A%TEB(3,l,m) = g1*A31 + g2*A32 + g3*A33
      A%Phi(1,l,m) = g1*A41 + g2*A42 + g3*A43 + g4*A44
    end do
  end do

end subroutine HealpixAlm_Sim_TEB


subroutine HealpixAlm_Sim(A, P, seed, HasPhi, dopol)
  use random
  use alm_tools
  use ran_tools
  Type(HealpixAlm) :: A
  Type(HealpixPower) :: P
  integer, intent(in), optional :: seed
  logical, intent(in), optional :: HasPhi,dopol
  integer l,m
  logical wantphi
  integer wantpol
  !namikawa
  !real(sp) xamp, corr, tamp, Bamp, Examp, sqrt2
  !complex(sp) g
  real(dp) xamp, corr, tamp, Bamp, Examp, sqrt2
  complex(dp) g

  if (present(seed)) then
    call InitRandom(seed)
  else
    if (.not. RandInited) call InitRandom
    RandInited = .true.
  end if

  wantpol = 1
  if (present(dopol)) then
    if (dopol) wantpol = 3
  end if

  if (present(HasPhi)) then
    wantphi= HasPhi
    if (wantphi .and. .not. associated(P%PhiCl)) call MpiStop('HealpixAlm_Sim: PhiCl not present')
  else
    wantphi = .false.
  end if

  call HealpixAlm_Init(A,P%lmax, wantpol, HasPhi=wantphi)
  sqrt2 = sqrt(2.)
  A%TEB = 0
  do l=1, P%lmax
    A%TEB(1,l,0) = Gaussian1()* sqrt(P%Cl(l,1))
    tamp = sqrt(P%Cl(l,1)/2)
    do m = 1, l
      A%TEB(1,l,m) =cmplx(Gaussian1(),Gaussian1())*tamp
    end do 
  end do

  if (wantphi) A%Phi=0
  if (wantpol >= 3) then  
    !polarization, E correlated to T
    do l = 2, P%lmax 
      if (p%cl(l,C_T) == 0) then
        tamp = 1.0  !Prevent divide by zero - TE should also be zero
      else
        tamp = p%cl(l,C_T)
      end if
      corr = p%cl(l,C_C)/tamp
      xamp = sqrt(P%cl(l,C_E) - corr*P%cl(l,C_C))
      Bamp = sqrt(P%cl(l,C_B))
      g = Gaussian1()
      if (wantphi) A%Phi(1,l,0) = g 
      A%TEB(2,l,0) = Corr*A%TEB(1,l,0) + real(g)*xamp
      A%TEB(3,l,0)=  Bamp*Gaussian1()
      xamp = xamp / sqrt2
      Bamp = Bamp /sqrt2
      do m =1, l
        g = cmplx(Gaussian1(),Gaussian1())
        A%TEB(2,l,m) = corr*A%TEB(1,l,m) + g*xamp
        A%TEB(3,l,m) = Bamp*cmplx(Gaussian1(),Gaussian1())     
        if (wantphi) A%Phi(1,l,m) = g 
      end do
    end do
  end if 

  if (wantphi) then
    !Make phi with correct correlation to T and E, AL May 2010
    do l=1, P%lmax
      if (P%Cl(l,1)==0) then
        tamp = 1.0
      else
        tamp=P%Cl(l,1)
      end if
      corr = P%PhiCl(l,2)/tamp
      if (wantpol >=3 .and. l>=2) then
        Examp = (P%PhiCl(l,3)-corr*P%cl(l,C_C))*sqrt( tamp/(p%cl(l,C_E)*tamp - p%cl(l,C_C)**2))
        xamp = sqrt(max(0._sp, P%PhiCl(l,1) - corr*P%PhiCl(l,2) - Examp**2 ))
        !write(*,*) l, tamp, p%cl(l,C_E)*tamp,p%cl(l,C_C)**2
        A%Phi(1,l,0) =  Examp * A%Phi(1,l,0)
        Examp = Examp/sqrt2
      else
        xamp = sqrt(max(0._sp,P%PhiCl(l,1) - corr*P%PhiCl(l,2)))
      end if        
      A%Phi(1,l,0) = A%Phi(1,l,0) + corr*A%TEB(1,l,0) + Gaussian1()*xamp
      xamp=  xamp/sqrt2
      do m = 1, l
        if (wantpol >=3 .and. l>=2) A%Phi(1,l,m) =  Examp * A%Phi(1,l,m) 
        A%Phi(1,l,m) = A%Phi(1,l,m) + corr*A%TEB(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp
      end do
    end do 
  end if

end subroutine HealpixAlm_Sim



  subroutine HealpixAlm_SimPhi(A, P, seed)
   use random
   use alm_tools
   use ran_tools
   Type(HealpixAlm) :: A
   Type(HealpixPower) :: P
   integer, intent(in), optional :: seed
   
   integer l,m

   if (present(seed)) then
      call InitRandom(seed)
   else
     if (.not. RandInited) call InitRandom
     RandInited = .true.
   end if
   call HealpixAlm_Init(A,P%lmax, 0,HasPhi = .true.)
   if (.not. P%lens) call MpiStop('must have phi power spectrum')

   A%Phi(:,0,:)=0
   do l=1, P%lmax
      A%Phi(1,l,0) =Gaussian1()* sqrt(P%PhiCl(l,1))
      do m = 1, l
       A%Phi(1,l,m) =cmplx(Gaussian1(),Gaussian1())* sqrt(P%PhiCl(l,1)/2)
      end do 
   end do

  end subroutine HealpixAlm_SimPhi


  subroutine HealpixMap_Init(M, npix, nmaps, nested, spinmap, HasPhi, pol, nside)

    Type(HealpixMap) :: M
    integer, intent(in), optional :: npix, nside 
    integer, intent(in), optional :: nmaps, spinmap
    logical, intent(in), optional :: nested, HasPhi, pol
    integer status

    if (present(npix)) then
     M%npix = npix
     if (present(nside)) then
      if (M%npix /= 12*nside**2) call MpiStop('HealpixMap_Init: nside and npix specified')
     end if
    else
     if (present(nside)) then
      M%npix = 12*nside**2
     else
      call MpiStop('HealpixMap_Init: must specifc nside or npix')
     end if  
    end if 

     call HealpixMap_Free(M)
        
     if (present(nmaps)) then
       M%nmaps = nmaps
     else
       M%nmaps = 1
     end if
     
     if (present(pol)) then
        if (pol) M%nmaps=3
       if (present(nmaps)) then
         if (nmaps>1) call MpiStop('HealpixMap_Init: currently only one pol map allowed')
       end if
     end if

     if (present(HasPhi)) then
       M%HasPhi = HasPhi
     else
       M%HasPhi = .false.
     end if

     if (M%nmaps /= 0) call HealpixMap_AllocateTQU(M,M%nmaps) 
 
     M%ordering = ord_ring
     if (present(nested)) then
        if (nested) M%ordering = ord_nest
     end if
  
     M%nside = npix2nside(M%npix)
     if (present(spinmap)) then
      M%spin = spinmap
     else
      M%spin = nospinmap
     end if
     if (M%spin /= nospinmap) then
       if (M%spin<1 .or. M%spin > 3) call MpiStop('Spin must be 0<spin<4')
       ALLOCATE(M%SpinField(0:M%npix-1),stat = status)
       M%spin = spinmap
     end if
    
    if (M%HasPhi) call HealpixMap_AllocatePhi(M)

  end subroutine HealpixMap_Init


  subroutine HealpixMap_AllocatePhi(M)
    Type(HealpixMap) :: M
    integer status
 
      ALLOCATE(M%Phi(0:M%npix-1),stat = status)
      if (status /=0) call MpiStop('HealpixMap_AllocatePhi: allocate')
      M%HasPhi = .true.

  end  subroutine HealpixMap_AllocatePhi

  subroutine HealpixMap_AllocateTQU(M, nmaps)
    Type(HealpixMap) :: M
    integer, intent(in) :: nmaps
    integer status
 
      M%nmaps = nmaps
      deallocate(M%TQU, stat=status)
      ALLOCATE(M%TQU(0:M%npix-1,M%nmaps),stat = status)
      if (status /= 0) call MpiStop('HealpixMap_AllocateTQU: allocate')
      M%TQU= 0
 
  end  subroutine HealpixMap_AllocateTQU


  subroutine HealpixMap_DeAllocateTQU(M)
    Type(HealpixMap) :: M
    integer status
 
    deallocate(M%TQU, stat=status)
    nullify(M%TQU)
    M%nmaps = 0

  end subroutine HealpixMap_DeAllocateTQU

  subroutine HealpixMap_PhiOnly(M)
    Type(HealpixMap) :: M
    integer status
       
   deallocate(M%TQU, stat=status)
   nullify(M%TQU)
   deallocate(M%SpinField, stat=status)
   nullify(M%SpinField)
   M%spin = nospinmap
   M%nmaps = 0   

  end subroutine HealpixMap_PhiOnly
  

  subroutine HealpixMap_Assign(MOut, Min)
   Type(HealpixMap) :: MOut, Min
   integer status

   call HealpixMap_Free(MOut)
   Mout = Min 
   nullify(MOut%TQU)
   if (Min%nmaps>0) then
     allocate(MOut%TQU(0:Min%npix-1,Min%nmaps),stat = status)
     if (status /= 0) call MpiStop('No Mem: HealpixMap_Assign')
     MOut%TQU = Min%TQU
   end if
   nullify(MOut%SpinField,MOut%phi)
   if (MIn%spin /= nospinmap) then
     allocate(MOut%SpinField(0:Min%npix-1),stat = status)
     if (status /= 0) call MpiStop('No Mem: HealpixMap_Assign')
     MOut%SpinField = Min%SpinField
     !MOut%SpinQ => MOut%SpinField(:,1)
     !MOut%SpinU => MOut%SpinField(:,2) 
   end if
   if (MIn%HasPhi) then
     allocate(MOut%Phi(0:Min%npix-1),stat = status)
     if (status /= 0) call MpiStop('No Mem: HealpixMap_Assign')
     MOut%Phi = Min%Phi   
   end if 

  end subroutine HealpixMap_Assign


  subroutine HealpixMap_Read(OutMAP,fname, map_limit)
   CHARACTER(LEN=80), DIMENSION(1:120) :: header_in
   character(LEN=*), intent(in) :: fname
   integer, optional :: map_limit
   integer nmaps
   
   Type(HealpixMap) OutMap, TmpMap 

   call HealpixMap_Free(OutMap)
    
   if (.not. FileExists(fname)) call MpiStop('HealpixMap_Read: File not found - '//trim(fname)) 
   OutMap%npix = getsize_fits(fname, nmaps=OutMap%nmaps, ordering=OutMap%ordering, nside=OutMap%nside,&
        type=OutMap%type)
   nmaps = outMap%nmaps
   if (present(map_limit)) then
   nmaps = min(nmaps, map_limit)
   end if
   if ((OutMap%ordering /=  ord_ring).and.(OutMap%ordering /= ord_nest)) then
     PRINT*,'The ordering scheme of the map must be RING or NESTED.'
     PRINT*,'No ordering specification is given in the FITS-header!'
     call MpiStop('')
    endif
  if (OutMap%nside /= npix2nside(OutMap%npix)) then ! accept only full sky map
     print*,'FITS header keyword NSIDE = ',OutMap%nside,' does not correspond'
     print*,'to the size of the map!'
     call MpiStop('')
   endif

   OutMap%HasPhi = .false.
   OutMap%spin = nospinmap
   if (OutMap%nmaps/=nmaps) then
     TmpMap=OutMap
     call HealpixMap_AllocateTQU(TmpMap,TmpMap%nmaps) 
     call input_map(fname, TmpMAP%TQU, OutMap%npix, TmpMap%nmaps, &
       &   fmissval=fmissval, header= header_in)
     OutMap%nmaps=nmaps
     nullify(OutMap%TQU)
     call HealpixMap_AllocateTQU(OutMap,nmaps)
     OutMap%TQU(:,1:nmaps)=TmpMap%TQU(:,1:nmaps)
     call HealpixMap_Free(TmpMap)     
   else
    call HealpixMap_AllocateTQU(OutMap,nmaps) 
    call input_map(fname, OutMAP%TQU, OutMap%npix, OutMap%nmaps, &
       &   fmissval=fmissval, header= header_in)
   end if
       
   
!!To do, boring...
  ! do j=1,nmaps
  !   call get_card(header_in, trim(numcat('TTYPE',j)), ttype(j))
  !   call get_card(header_in, trim(numcat('TUNIT',j)), tunit(j))
  ! enddo

  end subroutine HealpixMap_Read

  subroutine HealpixMap_Write(M, fname, overwrite)
    Type(HealpixMap), intent(in) :: M
    character(LEN=*), intent(in) :: fname
    logical, intent(in), optional :: overwrite
    CHARACTER(LEN=80), DIMENSION(1:120) :: header
    integer nlheader

  if (present(overwrite)) then
    if (overwrite) call DeleteFile(fname)
  else
    if (FileExists(fname)) call MpiStop('HealpixMap_Write: file already exists - '//trim(fname))
  end if

  header = ' '
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'COMMENT','     Sky Map Pixelisation Specific Keywords    ')
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'PIXTYPE','HEALPIX','HEALPIX Pixelisation')
  if (M%ordering == ord_ring) then
    call add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED')
  else
    call add_card(header,'ORDERING','NESTED',  'Pixel ordering scheme, either RING or NESTED')
  end if
  call add_card(header,'NSIDE'   ,M%nside,   'Resolution parameter for HEALPIX')
  call add_card(header,'FIRSTPIX',0,'First pixel # (0 based)')
  call add_card(header,'LASTPIX',M%npix-1,'Last pixel # (0 based)')
  call add_card(header) ! blank line
  call add_card(header,'CREATOR','HEALPixObj',        'Software creating the FITS file')
  call add_card(header,'INDXSCHM','IMPLICIT',' Indexing : IMPLICIT or EXPLICIT')
  call add_card(header,'GRAIN', 0, ' Grain of pixel indexing') ! full sky
  if (M%nmaps == 3) then
     call add_card(header,'POLAR',.true.," Polarisation included (True/False)")
  else
    call add_card(header,'POLAR',.false.," Polarisation included (True/False)")
  end if
  call add_card(header) ! blank line
  call add_card(header,"TTYPE1", "TEMPERATURE","Temperature map")
  call add_card(header,"TUNIT1", "muK", "map unit")
  call add_card(header)
  if (M%nmaps == 3) then
    call add_card(header,"TTYPE2", "Q-POLARISATION","Q Polarisation map")
    call add_card(header,"TUNIT2", "muK", "map unit")
    call add_card(header)
    call add_card(header,"TTYPE3", "U-POLARISATION","U Polarisation map")
    call add_card(header,"TUNIT3", "muK", "map unit")
    call add_card(header)
  end if
  call add_card(header,"COMMENT","*************************************")
  nlheader = SIZE(header)
  call write_bintab(M%TQU, M%npix, M%nmaps, header, nlheader, fname)
end subroutine HealpixMap_Write

  subroutine HealpixAlm_Write(A, fname)
    !Thanks to Sam Leach
    Type(HealpixAlm), intent(in) :: A
    character(LEN=*), intent(in) :: fname
    CHARACTER(LEN=80), DIMENSION(1:120) :: header
    integer nlheader,ii

    do ii=1,A%npol
       header = ' '
       call add_card(header,'COMMENT','-----------------------------------------------')
       call add_card(header,'COMMENT','     Sky Map Pixelisation Specific Keywords    ')
       call add_card(header,'COMMENT','-----------------------------------------------')
       if (ii == 1) then
          call add_card(header,"EXTNAME","""ANALYSED a_lms (TEMPERATURE)""")
       elseif (ii == 2) then
          call add_card(header,"EXTNAME","'ANALYSED a_lms (ELECTRIC component)'")
       elseif (ii == 3) then
          call add_card(header,"EXTNAME","'ANALYSED a_lms (CURL / MAGNETIC component)'")
       endif
       call add_card(header) ! blank line
       call add_card(header,'CREATOR','HEALPixObj',        'Software creating the FITS file')
       if (A%npol == 3) then
          call add_card(header,'POLAR',.true.," Polarisation included (True/False)")
       else
          call add_card(header,'POLAR',.false.," Polarisation included (True/False)")
       endif
       call add_card(header) ! blank line
       call add_card(header,"TTYPE1", "INDEX"," i = l^2 + l + m + 1")
       call add_card(header,"TUNIT1", "   "," index")
       call add_card(header)
       call add_card(header,"TTYPE2", "REAL"," REAL a_lm")
       call add_card(header,"TUNIT2", "muK"," alm units")
       call add_card(header)
        !
       call add_card(header,"TTYPE3", "IMAG"," IMAGINARY a_lm")
       call add_card(header,"TUNIT3", "muK"," alm units")
       call add_card(header)
       call add_card(header,"COMMENT","*************************************")
       
       nlheader = SIZE(header)
       call dump_alms(fname,A%TEB(ii,0:A%lmax,0:A%lmax),A%lmax,header,nlheader,ii-1)

    end do
  end  subroutine HealpixAlm_Write



  subroutine HealpixMap_Free(M)
   Type(HealpixMap) :: M
   integer status  

   call HealpixMap_DeAllocateTQU(M)

   !To prevent memory leaks, have to do fail-safe deallocates separately
   !otherwise can fail on the first one. Thanks Duncan Hanson!
   deallocate(M%SpinField, stat=status)
   deallocate(M%Phi, stat=status)
  
   nullify(M%SpinField,M%Phi)

  end subroutine HealpixMap_Free

  subroutine HealpixMap_Nullify(M)
   Type(HealpixMap) :: M
  
   nullify(M%TQU)
   nullify(M%SpinField,M%Phi)

  end subroutine HealpixMap_Nullify


  subroutine HealpixMap_ForceRing(M)
    USE pix_tools, ONLY : convert_nest2ring
   Type(HealpixMap) :: M
   integer i

   if (M%ordering /= ord_ring) then
      call convert_nest2ring (M%nside, M%TQU)
    ! do i=1, M%nmaps
    !  call convert_nest2ring (M%nside, M%TQU(:,i))
    ! end do
     if (M%spin/= nospinmap) then
       call MpiStop('ring not done for pol')
!       call convert_nest2ring (M%nside, M%SpinField(:,1))
!       call convert_nest2ring (M%nside, M%SpinField(:,2))
      end if
      if (M%HasPhi) then
       call convert_nest2ring (M%nside, M%Phi)
      end if
   end if
   M%ordering = ord_ring

  end subroutine HealpixMap_ForceRing

  
  subroutine HealpixMap_ForceNest(M)
    USE pix_tools, ONLY : convert_ring2nest
   Type(HealpixMap) :: M
   integer i

   if (M%ordering /= ord_nest) then
     call convert_ring2nest (M%nside, M%TQU)
    ! do i=1, M%nmaps
    !  call convert_ring2nest (M%nside, M%TQU(:,i))
    ! end do
     if (M%spin/= nospinmap) then
      call MpiStop('nest not done for pol')
       !call convert_ring2nest (M%nside, M%SpinField(:,1))
       !call convert_ring2nest (M%nside, M%SpinField(:,2))
      end if

      if (M%HasPhi) then
       call convert_ring2nest (M%nside, M%Phi)
      end if
     end if
   M%ordering = ord_nest

  end subroutine HealpixMap_ForceNest


  subroutine HealpixMapMulCut(InMap,CutMap,OutMap, map_ix, missval)
    Type(HealpixMap), intent(in) :: InMap,CutMap
    Type(HealpixMap), intent(out) :: OutMap
    integer, intent(in), optional :: map_ix
    real(sp), intent(in), optional :: missval
    integer i, j, ix

   if (present(map_ix)) then
     ix = map_ix
   else
     ix = 2
   end if
   call HealpixMap_ForceRing(InMap)
   call HealpixMap_ForceRing(CutMap)
   if (InMap%npix /= CutMap%npix) call MpiStop('HealpixMapMulCut: Map size mismatch')
   call HealpixMap_Init(OutMap,InMap%npix,InMap%nmaps)
   outMap%ordering  = ord_ring
   if (CutMap%nmaps < ix ) call MpiStop('HealpixMapMulCut: not enough maps')
   if (present(missval)) then
 
      do j=0, InMap%npix -1
        if (CutMap%TQU(j,ix) == 0) then
          OutMap%TQU(j,:) = missval
        else       
          OutMap%TQU(j,:) = InMap%TQU(j,:)
        end if
     end do
 
   else

    do i=1,InMap%nmaps
     do j=0, InMap%npix -1
     OutMap%TQU(j,i) = InMap%TQU(j,i) * CutMap%TQU(j,ix)
     end do
    end do

   end if
  end subroutine HealpixMapMulCut

   subroutine HealpixMap_MulCutFile(InMap,CutFile,OutMap)
    Type(HealpixMap), intent(in) :: InMap
    Type(HealpixMap), intent(out) :: OutMap
    character(len=*), intent(in) :: CutFile
    Type(HealpixMap) CutMap

    call HealpixMap_Read(CutMap,Cutfile)
    call HealpixMapMulCut(InMap,CutMap,OutMap)
    call HealpixMap_Free(CutMap)

   end subroutine HealpixMap_MulCutFile


   subroutine HealpixMap2alm(H, M,A, almax,theta_cut_deg,map_ix, dopol)

     Type (HealpixInfo) :: H
     Type(HealpixMap), intent(in) :: M
     integer, intent(in) :: almax
     real(dp), intent(in), optional :: theta_cut_deg
      integer, intent(in), optional :: map_ix
      logical, intent(in), optional :: dopol
     Type(HealpixAlm) :: A
     integer npol, ix
     real(dp) cos_theta_cut

     call HealpixMap_ForceRing(M)

     npol = 1
     if (present(dopol)) then
      if (dopol) npol =3
     else
      if (M%nmaps ==3) npol = 3
     end if

     if (M%nmaps ==0) npol = 0

     if (present(theta_cut_deg)) then
       cos_theta_cut =  SIN(theta_cut_deg/180.d0*HO_pi)
       if (theta_cut_deg < 0) cos_theta_cut = -1
      else
       cos_theta_cut = -1
     end if

     call HealpixAlm_Init(A,almax, npol, M%spin, M%HasPhi)
 
     if (npol==1 .and. M%nmaps>0) then
      ix = 1
      if (present(map_ix)) ix = map_ix
      call map2scalalm(H, almax, M%TQU(:,ix), A%TEB,cos_theta_cut)
     else if (npol ==3 .and. M%nmaps>0) then
      if (present(map_ix)) call MpiStop(' cannot have polarization and multiple map indices')
      call map2polalm(H, almax, M%TQU, A%TEB,cos_theta_cut)
     end if  

     if (M%HasPhi) then
       call map2scalalm(H, almax, M%Phi, A%Phi,cos_theta_cut)
     end if

     if (M%spin /= nospinmap) then
       !Haven't computed ring weights for odd spins
        call map2spinalm(H,almax, M%SpinField,  A%SpinEB, M%spin,cos_theta_cut)
     end if

end subroutine HealpixMap2alm

subroutine HealpixMap2alm_nami(H, M,A, almax,theta_cut_deg,map_ix, dopol)

     Type (HealpixInfo) :: H
     Type(HealpixMap), intent(in) :: M
     integer, intent(in) :: almax
     real(dp), intent(in), optional :: theta_cut_deg
      integer, intent(in), optional :: map_ix
      logical, intent(in), optional :: dopol
     Type(HealpixAlm) :: A
     integer npol, ix
     real(dp) cos_theta_cut

     call HealpixMap_ForceRing(M)

     npol = 1
     if (present(dopol)) then
      if (dopol) npol =3
     else
      if (M%nmaps ==3) npol = 3
     end if

     if (M%nmaps ==0) npol = 0

     if (present(theta_cut_deg)) then
       cos_theta_cut =  SIN(theta_cut_deg/180.d0*HO_pi)
       if (theta_cut_deg < 0) cos_theta_cut = -1
      else
       cos_theta_cut = -1
     end if

     call HealpixAlm_Init(A,almax, npol, M%spin, M%HasPhi)
 
      if (present(map_ix)) call MpiStop(' cannot have polarization and multiple map indices')

 
end subroutine HealpixMap2alm_nami

   subroutine HealpixMapArray_Free(M)
     Type(HealpixMap) :: M(:)
     integer i
     
     do i=1, size(M) 
      call HealpixMap_Free(M(i))
     end do
     
   end subroutine HealpixMapArray_Free

   subroutine HealpixMapSet2CrossPowers(H, M, Pows, nmap, almax, dofree)
     integer, intent(in) :: nmap
     Type (HealpixInfo) :: H
     Type(HealpixMap), intent(in) :: M(nmap)
     Type(HealpixCrossPowers) :: Pows
     integer, intent(in) :: almax
     integer i
     logical :: dofree
     Type(HealpixMapArray) :: maps(nmap)

!Does not deallocate pows, assumed undefined

     if (nmap<0)  call MpiStop('HealpixMapSet2CrossPowers: must have one or more maps')
     do i=1,nmap
      call HealpixMap_ForceRing(M(i))
     end do
     
     Pows%nmaps = nmap
     Pows%lmax = almax
     Pows%npol=0
     do i=1,nmap
       if (M(i)%nmaps /=1 .and. M(i)%nmaps /=3) call MpiStop('HealpixMapSet2CrossPowers: must be scalar or pol')
        Pows%npol = max(Pows%npol,M(i)%nmaps)
        if(dofree) then
         allocate(maps(i)%M(0:size(M(i)%TQU,1)-1,size(M(i)%TQU,2)))
         maps(i)%M = M(i)%TQU
         call HealpixMap_Free(M(i))
        else
        maps(i)%M => M(i)%TQU
        end if
     end do
 
 
     if (Pows%npol==1) then
      call maparray2scalcrosspowers(H, almax, maps, Pows,  nmap, dofree)
     else
      call maparray2crosspowers(H, almax, maps, Pows,  nmap, dofree)
     end if  
 
   end subroutine HealpixMapSet2CrossPowers


  subroutine HealpixMap_Smooth(H, MapIn, MapOut, lmax, fwhm)
       Type(HealpixInfo) :: H
       Type(HealpixMap) :: MapIn, MapOut
       integer, intent (in) :: lmax
       Type(HealpixAlm) :: A
       real(dp), intent(in) :: fwhm
       
      call HealpixMap2Alm(H,MapIn, A, lmax)
      call HealpixAlm_Smooth(A, fwhm)       
      call HealpixAlm2Map(H,A, MapOut, MapIn%npix)
      call HealPixAlm_Free(A)
      
  end subroutine HealpixMap_Smooth  

  subroutine HealpixMap_SimulateUnlensed(H, M, P, Beam, want_pol)
       Type(HealpixInfo) :: H
       Type(HealpixMap) :: M
       Type(HealpixPower) :: P, BeamP
       real(dp), intent(in) :: Beam(0:)
       logical, intent(in) :: want_pol
       Type(HealpixAlm) :: A
       
       call HealpixPower_Assign(BeamP,P)
       call HealpixPower_Smooth_Beam(BeamP, Beam,-1)
       call HealpixAlm_Sim(A, BeamP, HasPhi=.false., dopol = want_pol)
       call HealpixAlm2Map(H, A, M, nside2npix(H%nside))
       call HealpixAlm_Free(A)
       call HealpixPower_Free(BeamP)
      
  end subroutine HealpixMap_SimulateUnlensed


  subroutine HealpixAlm2GradientMap(H, A, M, npix, What)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: M
     Type(HealpixAlm), intent(in) :: A
     Type(HealpixAlm) :: AT
     integer, intent(in) :: npix
     character(LEN=*), intent(in) :: What
   
    call HealpixMap_Init(M,npix,nmaps = 0, spinmap = 1)
    if (What(1:1) == 'P') then
     if (.not. A%HasPhi) call MpiStop('HealpixAlm2GradientMap: No phi field')
     call alm2GradientMap(H, A%lmax, A%Phi,M%SpinField)
    else if (What(1:1) == 'T') then
     call HealpixAlm_Init(AT, A%lmax,npol = 0, HasPhi = .true.)
     AT%Phi = A%TEB(1:1,:,:)
     call alm2GradientMap(H,A%lmax, AT%Phi,M%SpinField)
     call HealpixAlm_Free(AT)
    else
     call MpiStop('HealpixAlm2GradientMap: unknown field')
    end if
end subroutine HealpixAlm2GradientMap


subroutine  HealpixExactLensedMap(H,A, M, npix)
  Type (HealpixInfo) :: H
  Type(HealpixMap) :: GradPhi, M
  integer, intent(in) :: npix
  Type(HealpixAlm), intent(in) :: A

  write(*,*) "alm to gradient map"
  call HealpixAlm2GradientMap(H,A,GradPhi,npix,'PHI')
  write(*,*) "compute gradient phi"
  call HealpixExactLensedMap_GradPhi(H,A, GradPhi,M)
  write(*,*) "free memory"
  call HealpixMap_Free(GradPhi)

end subroutine  HealpixExactLensedMap


subroutine HealpixExactLensedMap_GradPhi(H,A, GradPhi,M)
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: GradPhi, M
  Type(HealpixAlm), intent(in) :: A

  call HealpixMap_Init(M,GradPhi%npix,nmaps = A%npol)
  if (GradPhi%spin /=1) call MpiStop('HealpixExactLensedMap: GradPhi must be spin 1 field')
  if (A%npol ==1) then
    write(*,*) "alm to lensed map"
    call scalalm2LensedMap(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU(:,1))
  else
    call alm2LensedMap(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU)
  end if
end subroutine HealpixExactLensedMap_GradPhi


subroutine HealpixInterpLensedMap(H,A, M, npix, factor, interp_method)
  implicit none
  Type (HealpixInfo) :: H
  Type(HealpixMap) :: GradPhi, M
  integer, intent(in) :: npix
  Type(HealpixAlm), intent(in) :: A
  real, intent(in), optional :: factor
  integer, intent(in), optional :: interp_method
  real fact
  integer method
     
  fact = 1.5 
  if (present(factor)) fact = factor
  method = interp_basic
  if (present(interp_method)) method = interp_method
  write(*,*) "alm to map (gradient)"
  call HealpixAlm2GradientMap(H,A,GradPhi,npix,'PHI')
  write(*,*) "lensed map"
  call HealpixInterpLensedMap_GradPhi(H,A,GradPhi,M,fact,interp_method)
  write(*,*) "free memory"
  call HealpixMap_Free(GradPhi)

end subroutine HealpixInterpLensedMap


subroutine HealpixInterpLensedMap_GradPhi(H,A,GradPhi,M,factor,interp_method)
  implicit none
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: GradPhi, M
  Type(HealpixAlm), intent(in) :: A
  real, intent(in) :: factor
  integer, intent(in) :: interp_method
      
  call HealpixMap_Init(M,GradPhi%npix,nmaps = A%npol)
  if(GradPhi%spin/=1) call MpiStop('HealpixExactLensedMap: GradPhi must be spin 1 field')
  if(interp_method==interp_basic) then
    if(abs(factor-nint(factor))>1e-4) call MpiStop('interp_factor must be 2^n for interp_basic')
  end if
  if(A%npol ==1) then
    if(interp_method==interp_basic) then
      call scalalm2LensedmapInterp(H,A%lmax,A%TEB,GradPhi%SpinField,M%TQU,nint(factor))
    else
      call scalalm2LensedmapInterpCyl(H,A%lmax,A%TEB,GradPhi%SpinField,M%TQU, factor)
    end if 
  else
    if(interp_method==interp_basic) then
      call alm2LensedmapInterp(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU, nint(factor))
    else
      write(*,*) '!'
      call alm2LensedmapInterpCyl(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU, factor)
    end if
  end if

end subroutine HealpixInterpLensedMap_GradPhi


  subroutine  HealpixQuadLensedMap(H,A, M, npix)
   ! ** Note does not give sky with accurate lensed C_l at l>~1200 **
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     integer, intent(in) :: npix
     Type(HealpixAlm), intent(in) :: A
      
     call HealpixAlm2GradientMap(H,A,GradPhi,npix,'PHI')
     call HealpixQuadLensedMap_GradPhi(H,A, GradPhi,M)
     call HealpixMap_Free(GradPhi)

  end subroutine  HealpixQuadLensedMap


  subroutine  HealpixQuadLensedMap_GradPhi(H,A, GradPhi,M)
   ! ** Note does not give sky with accurate lensed C_l  at l>~1200 **
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     Type(HealpixAlm), intent(in) :: A
      
     call HealpixMap_Init(M,GradPhi%npix,nmaps = A%npol)
     if (GradPhi%spin /=1) call MpiStop('HealpixExactLensedMap: GradPhi must be spin 1 field')
     if (A%npol ==1) then
      call alm2LensedQuadContrib(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU(:,1))
     else
     call MpiStop('not done yet')
     end if
  end subroutine  HealpixQuadLensedMap_GradPhi
   

  subroutine HealpixAlm2Map(H, A, M, npix, DoPhi, DoT)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: M
     Type(HealpixAlm), intent(in) :: A
     integer, intent(in) :: npix
     logical, intent(in), optional :: DoPhi, DoT
     logical Phi
     integer npol

     if (present(DoPhi)) then
      Phi = DoPhi .and. A%HasPhi
     else
      Phi = A%HasPhi
      end if

     if (present(DoT)) then
       if (DoT) then 
           npol = A%npol 
       else  
           npol = 0 
       end if
     else
       npol = A%npol
      end if
       
    call HealpixMap_Init(M,npix,npol, spinmap = A%spin, HasPhi = Phi)

    if (M%nmaps > 0) then
     if (npol>1) then
        call polalm2map(H,A%lmax, A%TEB,M%TQU)
      else
        call scalalm2map(H,A%lmax, A%TEB,M%TQU(:,1))
      end if    
    end if
    if (M%spin /= nospinmap) then
        call spinalm2map(H,A%lmax, A%SpinEB,M%SpinField, A%spin)
    end if
    if (Phi) then
        call scalalm2map(H,A%lmax, A%Phi,M%Phi)
    end if

  end subroutine HealpixAlm2Map

  subroutine HealpixMap_PolToSpin2Field(M)
      Type(HealpixMap) :: M
      integer status

       deallocate(M%SpinField,stat =status)
       ALLOCATE(M%SpinField(0:M%npix-1),stat = status)
       M%spin = 2
       M%SpinField = cmplx(M%TQU(:,2),M%TQU(:,3))

  end  subroutine HealpixMap_PolToSpin2Field

  subroutine HealpixMap_Spin2FieldTopol(M, delspin)
      Type(HealpixMap) :: M
      logical, intent(in), optional :: delspin

       if (M%nmaps==0) then
         call HealpixMap_AllocateTQU(M,3) 
       end if
       M%TQU(:,2) = real(M%SpinField)
       M%TQU(:,3) = aimag(M%SpinField)
       if (present(delspin)) then
       if (delspin) then
          deallocate(M%SpinField)
           M%spin = nospinmap
        end if
       end if

  end  subroutine HealpixMap_Spin2FieldToPol


 subroutine HealpixMap_SetToIndexOnly(M, ix)
    Type(HealpixMap) :: M
    integer, intent(in) :: ix
    REAL(SP), DIMENSION(:,:), allocatable :: TQU
    
    if (M%nmaps < ix) call MpiStop('HealpixMap_SetToIndexOnly: index out of bounds')
    M%nmaps =1
    allocate(TQU(0:M%npix-1,1))
    TQU(:,1) = M%TQU(:,ix)
    call HealpixMap_AllocateTQU(M,1) 
    M%TQU = TQU
    deallocate(TQU)
    
 end subroutine HealpixMap_SetToIndexOnly
 
  subroutine HealpixMap_SetToDataOnly(M, want_pol)
    Type(HealpixMap) :: M
    REAL(SP), DIMENSION(:,:), allocatable :: TQU
    logical, intent(in) :: want_pol
    
    if (.not. want_pol) then
     call HealpixMap_SetToIndexOnly(M,1)
    else 
    if (M%nmaps <3) call MpiStop('HealpixMap_SetToDataOnly: not enough maps')
    M%nmaps =3
    allocate(TQU(0:M%npix-1,3))
    TQU(:,1:3) = M%TQU(:,1:3)
    call HealpixMap_AllocateTQU(M,3) 
    M%TQU = TQU
    deallocate(TQU)
    end if
    
 end subroutine HealpixMap_SetToDataOnly
 
  subroutine HealpixMap_AddPol(M)
    Type(HealpixMap) :: M
    REAL(SP), DIMENSION(:), allocatable :: T
    
    if (M%nmaps >1 ) call MpiStop('HealpixMap_AddPol:already more than one map')
    M%nmaps =3
    allocate(T(0:M%npix-1))
    T = M%TQU(:,1)
    call HealpixMap_AllocateTQU(M,3) 
    M%TQU(:,1) = T
    deallocate(T)
    
 end subroutine HealpixMap_AddPol
 
  subroutine HealpixMap_udgrade(M, Mout, nside_out, pessimistic)
     use udgrade_nr
     Type(HealpixMap) :: M, Mout
     integer, intent(in) :: nside_out
     logical, intent(in), optional :: pessimistic
     integer i
     logical isring, pess
     

    if (present(pessimistic)) then 
      pess = pessimistic
    else
      pess = .false.
    end if
    if (nside_out== M%nside) then
     call HealpixMap_Assign(Mout, M)
    else
    isring = M%ordering == ord_ring
    call HealpixMap_ForceNest(M)
    call HealpixMap_Init(Mout,nside2npix(nside_out), M%nmaps, nested=.true.)
     do i=1, M%nmaps
       call sub_udgrade_nest(M%TQU(:,i),M%nside,MOut%TQU(:,i),nside_out, fmissval, pess)
    end do
    if (isring) call HealpixMap_ForceRing(MOut)
    end if
  end subroutine HealpixMap_udgrade

  subroutine HealpixMap_HaarTransform(M,Mdegrade,Mdetail)
   implicit none
   Type(HealpixMap) :: M, Mdegrade, Mdetail
   integer i, j
   real :: bas1(4) = (/ -1, -1,  1,  1 /) /4.
  ! real :: bas2(4) = (/ -1,  1,  0,  0 /) /2.
  ! real :: bas3(4) = (/  0,  0, -1,  1 /) /2.
   real :: bas2(4) = (/ -1,  1,  1,  -1 /)  /4.
   real :: bas3(4) = (/ -1,  1, -1,   1 /) /4.
   if (M%npix <=12) call MpiStop( 'Map too coarse to Haar transform')
   call HealpixMap_ForceNest(M)
   call HealpixMap_Init(Mdetail,M%npix/4, M%nmaps*3, nested = .true.)
   call HealpixMap_udgrade(M,Mdegrade,M%nside/2, pessimistic = .true.)
   
   do j=1, M%nmaps
    do i=0,Mdetail%npix -1
      if (any(M%TQU(i*4:i*4+3,j) == fmissval)) then
        Mdetail%TQU(i,(j-1)*4+1:(j-1)*4+3) = fmissval
      else
       Mdetail%TQU(i,(j-1)*4+1) = sum(M%TQU(i*4:i*4+3,j)*bas1)
       Mdetail%TQU(i,(j-1)*4+2) = sum(M%TQU(i*4:i*4+3,j)*bas2)
       Mdetail%TQU(i,(j-1)*4+3) = sum(M%TQU(i*4:i*4+3,j)*bas3)
      end if
    end do
   end do

  end subroutine HealpixMap_HaarTransform


  subroutine HealpixMap_HaarReconstruct(Mdegrade, Mdetail, M)
     Type(HealpixMap) :: M, Mdegrade, Mdetail
   integer i, j
   real :: bas1(4) = (/ -1, -1,  1,  1 /) 
!   real :: bas2(4) = (/ -1,  1,  0,  0 /) 
!   real :: bas3(4) = (/  0,  0, -1,  1 /) 
   real :: bas2(4) = (/ -1,  1,  1,  -1 /) 
   real :: bas3(4) = (/ -1,  1, -1,  1 /) 

   if (Mdegrade%npix /= Mdetail%npix) call MpiStop('map size mismatch')
   if (mod(Mdetail%nmaps,3)/=0) call MpiStop('detail not sets of 3 maps') 
   call HealpixMap_ForceNest(Mdegrade)
   call HealpixMap_ForceNest(Mdetail)
   call HealpixMap_Init(M,Mdegrade%npix*4, Mdegrade%nmaps/3, nested = .true.)
   call HealpixMap_udgrade(Mdegrade,M,Mdegrade%nside*2)
   
   do j=1, M%nmaps
    do i=0,Mdetail%npix -1
        if (any(M%TQU(i*4:i*4+3,j) == fmissval)) then
          M%TQU(i*4:i*4+3,j) = fmissval
        else
         M%TQU(i*4:i*4+3,j) = M%TQU(i*4:i*4+3,j) + Mdetail%TQU(i,(j-1)*4+1)*bas1 &
          + Mdetail%TQU(i,(j-1)*4+2)*bas2 + Mdetail%TQU(i,(j-1)*4+3)*bas3
        end if
    end do
   end do

  end subroutine HealpixMap_HaarReconstruct

  subroutine HaarComponents_Free(C) 
    Type(HaarComponents) :: C
    integer i

    if (associated(C%details)) then
     do i=1, C%order
      call HealpixMap_Free(C%details(i))
     end do
     deallocate(C%details)
     nullify(C%details)
    end if
    call HealpixMap_Free(C%degraded)
 
  end subroutine HaarComponents_Free

  subroutine HealpixMap_HaarComponents(M,C, order)
     Type(HealpixMap) :: M, D
     Type(HaarComponents) :: C
     integer, intent(in) :: order
     integer i
  
     C%order = order
     allocate(C%details(order))
     call HealpixMap_Assign(D,M)
     do i=1, order
       call HealpixMap_HaarTransform(D,C%degraded,C%details(i))     
       if (i/= order) call HealpixMap_Assign(D, C%degraded)
     end do
  end subroutine HealpixMap_HaarComponents

  
   subroutine HealpixMap_FromHaarComponents(C,M)
     Type(HealpixMap) :: M, D
     Type(HaarComponents) :: C
     integer i
  
     call HealpixMap_Assign(D,C%Degraded)
     do i=C%order,1,-1
       call HealpixMap_HaarReconstruct(D, C%details(i), M)
       if (i/=1) call HealpixMap_Assign(D, M)
     end do

  end subroutine HealpixMap_FromHaarComponents


 subroutine HealpixMap_HarrPowerMap(C, M, nside_power, nside_map)
!Map the power in the three details making up each pixel is nside_power averaged at resolution 
!nside_map

   Type(HaarComponents) :: C
   Type(HealpixMap) :: M
   integer, intent(in) :: nside_power, nside_map
   integer i, fac, j, pixcount
  

   if (C%degraded%nmaps /= 1) call MpiStop('Only does power spectra for single map')
   if (nside_power < nside_map) call MpiStop( 'Can only make map of power at lower resolution')

   call HealpixMap_Init(M, nside2npix(nside_map), nested = .true.)
   fac = nside2npix(nside_power) / M%npix
   j = 1
   do while (C%details(j)%nside > nside_power)
    j=j+1
    if (j> C%order) call MpiStop('Haar coeffs not available')
   end do
   if (C%details(j)%nside /= nside_power) call MpiStop('Haar coeffs not available')
   do i=0, M%npix-1
     pixcount =  count( C%details(j)%TQU(i*fac:(i+1)*fac - 1,:) /= fmissval) 
     if (pixcount /= 0) then
      M%TQU(i,1) = sum(C%details(j)%TQU(i*fac:(i+1)*fac - 1,:)**2, &
              mask =  (C%details(j)%TQU(i*fac:(i+1)*fac - 1,:) /= fmissval) ) / pixcount
     else
      M%TQU(i,1) = fmissval
     end if         
   end do      

  
 end subroutine HealpixMap_HarrPowerMap

  subroutine HealpixMap_HarrPowerSpec(C, P)
   Type(HaarComponents) :: C
   real(dp) P(C%order)
   integer i

   if (C%degraded%nmaps /= 1) call MpiStop('Only does power spectra for single map')
   do i=1, C%order
    P(i) = sum(C%details(i)%TQU(:,:)**2, mask = C%details(i)%TQU(:,:) /= fmissval) / &
           count(C%details(i)%TQU(:,:)/=fmissval)
   end do      
  
 end subroutine HealpixMap_HarrPowerSpec

 subroutine CrossPowersToHealpixPowerArray(CrossPowers,PowerArray, dofree)
  Type(HealpixCrossPowers) :: CrossPowers
  Type(HealpixPower) :: PowerArray(*)
  logical, intent(in), optional :: dofree
  integer i,j,ix
  
    if (CrossPowers%npol>1) call MpiStop('CrossPowersToHealpixPowerArray: not done for pol')
    ix=0
      do i=1, CrossPowers%nmaps
        do j=1,i
        ix= ix+1
        call HealpixPower_Init(PowerArray(ix), CrossPowers%lmax, CrossPowers%npol==3, dolens=.false., nofree=.true.)
        PowerArray(ix)%Cl(:,C_T) = CrossPowers%Ps(i,j)%Cl(:,1,1)
        if (present(dofree)) then
         if (dofree) deallocate(CrossPowers%Ps(i,j)%Cl)
        end if        
       end do
    end do   
    if (dofree) deallocate(CrossPowers%Ps)

 end subroutine CrossPowersToHealpixPowerArray

end module HealpixObj
