module HealpixVis
 !Code largely courtesy of Mark Ashdown
 !Visualise the sky looking at sphere from the centre
  use HealpixObj  
  use Healpix_Types
 ! use AMLUtils
  implicit none

  ! visualisation PPM type
  logical :: VisDebugMsgs = .true.
  logical :: HealpixVis_force_range = .false.
  real :: HealpixVis_abs_max = 500.
  integer, parameter :: HealpixVis_DEF_WIDTH = 800, HealpixVis_DEF_HEIGHT=400
  
  type HealpixPPM
    integer :: nx, ny
    integer, dimension(:,:,:), pointer :: rgb  => NULL() 
  end type HealpixPPM

  type HealpixCMap
    integer :: n
    real, dimension(:), pointer :: x  => NULL() 
    real, dimension(:,:), pointer :: rgb  => NULL() 
  end type HealpixCMap

  real, dimension(3), parameter, private :: &
      ppm_background = (/ 1.0, 1.0, 1.0 /), &
      ppm_null = (/ 1.0, 1.0, 1.0 /)
   !  ppm_background = (/ 1.0, 0., 0. /)

  integer, parameter :: splot_none=0, splot_amp=1, splot_angle =2, splot_lines =3, &
        splot_Q =4, splot_U=5, plot_phi = -1

contains

  !======================================================================

  ! HealpixPPM constructor

  subroutine HealpixVis_ppm_init(ppm, nx, ny)
    type(HealpixPPM), intent(inout) :: ppm
    integer, intent(in) :: nx, ny
    integer err
    ppm%nx = nx
    ppm%ny = ny
    deallocate(ppm%rgb,stat = err)
    allocate(ppm%rgb(3, nx, ny))
  end subroutine HealpixVis_ppm_init

  !======================================================================

  ! HealpixPPM destructor

  subroutine HealpixVis_ppm_free(ppm)
    type(HealpixPPM), intent(inout) :: ppm
    ppm%nx = 0
    ppm%ny = 0
    deallocate(ppm%rgb)   
    nullify(ppm%rgb) 
  end subroutine HealpixVis_ppm_free

  !======================================================================

  ! Write PPM image to file

  subroutine HealpixVis_ppm_write(ppm, filename)

    type(HealpixPPM), intent(in) :: ppm
    character(len=*), intent(in) :: filename

    integer :: i, j, k

    open(1, file=filename, action='write', status='replace', &
        recl=3*ppm%nx*ppm%ny+200)
    write(1, '("P6")')
    write(1, '("# Created by HealpixVis")')
    write(1, '(i4," ",i4)') ppm%nx, ppm%ny
    write(1, '(i4," ")',advance='no') 255
  
    do k = 1 , ppm%ny
      do j = 1, ppm%nx
        do i = 1, 3
          write(1, '(a)',advance='no') char(ppm%rgb(i, j, k))
        end do
      end do
    end do
    close(1)

  end subroutine HealpixVis_ppm_write

  !======================================================================

  ! Plot map as PPM image using given projection and colourmap

  subroutine HealpixVis_map2ppm(M, ppm, projection, n, nx, ny, plot, symmetric, size_scale)
    type(HealpixMap), intent(in) :: M
    type(HealpixPPM), intent(inout) :: ppm
    integer, intent(in), optional :: n
    integer, intent(in), optional :: nx, ny
    integer, intent(in), optional :: plot
    logical, intent(in), optional :: symmetric
    real(SP), intent(in), optional :: size_scale
    character(LEN=*), intent(in), optional :: projection
    character proj
    integer plotopt
    integer nn, nyy,nxx
    logical symm
    real(SP) scale
    real(SP), dimension(:), allocatable :: AmpArr

    if (present(plot)) then
      plotopt = plot
     else
       plotopt = splot_none
    end if   
    if (present(size_scale)) then
     scale = size_scale
    else 
     scale=1 
    end if

    if (present(symmetric)) then
     symm = symmetric
    else 
     symm = .false. 
    end if
    if (present(Projection)) then
     proj = projection(1:1)
    else
     proj = 'M'
    end if


     if (present(n)) then
       nn = n
      else
       nn=1
      end if

     if (present(nx)) then
         nxx = nint(nx*scale)
        else
         if (proj=='O') then
           nxx = nint(HealpixVis_DEF_HEIGHT * scale)
          else
           nxx = nint(HealpixVis_DEF_WIDTH * scale)
          end if
       end if

       if (present(ny)) then
         nyy = nint(ny*scale)
        else
         nyy = nint(HealpixVis_DEF_HEIGHT*scale)
       end if

    if (plotopt /= splot_none) then
      if (M%Spin== 0 .and. plotopt /= plot_phi) stop 'HealpixVis_map2ppm: No spin field!'
        allocate(AmpArr(0:M%npix-1))
        if (plotopt == splot_amp) then
            AmpArr = abs(M%SpinField(:))
        else if (plotopt == splot_angle) then
            AmpArr = atan2(aimag(M%SpinField(:)),real(M%SpinField(:)))/M%spin
        else if (plotopt == splot_Q) then
            AmpArr = real(M%SpinField(:))
        else if (plotopt == splot_U) then
            AmpArr = aimag(M%SpinField(:))
        else if (plotopt == plot_phi) then
            AmpArr = M%Phi
        end if
        if (proj == 'M') then
          call HealpixVis_map2ppm2(M, ppm, HealpixVis_proj_mol,  nn, nxx, nyy, AmpArr, symmetric=symm)
        else if (proj=='O') then
         call HealpixVis_map2ppm2(M, ppm, HealpixVis_proj_orth,  nn, nxx, nyy, AmpArr, symmetric=symm)
        else
          stop 'HealpixVis_map2ppm: unsupported projection'
        end if
        deallocate(AmpArr)
    else
        if (proj == 'M') then
          call HealpixVis_map2ppm2(M, ppm, HealpixVis_proj_mol,  nn, nxx, nyy, symmetric=symm)
        else if (proj=='O') then
         call HealpixVis_map2ppm2(M, ppm, HealpixVis_proj_orth,  nn, nxx, nyy, symmetric=symm)
        else
          stop 'HealpixVis_map2ppm: unsupported projection'
        end if
    end if 

  end subroutine HealpixVis_map2ppm

  subroutine HealpixVis_Map2Python(M, fname, nsideview, map_ix)
      type(HealpixMap), intent(in) :: M
      integer, intent(in), optional :: nsideview 
      character (LEN=*), intent(in) :: fname
      integer, intent(in), optional :: map_ix
      type(HealpixMap) :: Mout 
      integer pix, nside
      real(dp) verts(3,4)
      real color(3), aval

      type(HealpixCMap) :: cmap
      real mapmin, maprange
      integer pol

      if (present(map_ix)) then
       pol = map_ix
      else      
       pol = 1 
      end if

      if (present(nsideview)) then
        nside = nsideview
      else
        nside = M%nside
      end if

     call HealpixMap_udgrade(M,Mout, nside)
     call HealpixVis_cmap_rainbow(cmap)
     mapmin = minval(Mout%TQU(:,pol))
     maprange = maxval(MOut%TQU(:,pol))-mapmin

     open(1, file=fname, action='write', status='replace')
     write (1,*) Mout%npix
     do pix = 0, Mout%npix-1 
       aval = Mout%TQU(pix,pol)
      ! if (aval /= fmissval) then
        if (aval==0 .or. (aval == fmissval) ) then
         color = 0
        else
         color = HealpixVis_getcol(cmap,(aval-mapmin)/maprange)
        end if
       call HealpixMap_Pix2Vertex(Mout,pix,verts)
       write (1,'(15f10.6)') verts(:,1),verts(:,2),verts(:,3),verts(:,4), color
     ! end if
     end do

     close(1)
     call HealpixMap_Free(Mout)
     call HealpixVis_cmap_free(cmap)
 

  end subroutine HealpixVis_Map2Python

  subroutine HealpixVis_map2Anim(M,fileroot,map_ix)
   type(HealpixMap), intent(in) :: M
   character(LEN=*), intent(in) :: fileroot
   integer, intent(in), optional :: map_ix
   type(HealpixPPM) :: ppm
   real(dp) :: R(3,3), ang
   integer pol
   integer, parameter :: nframes = 10
   character(LEN=3) :: frame
   integer ix

   stop 'this doesn''t work'
   if (present(map_ix)) then
    pol = map_ix
   else
    pol =1
    end if
   do ix = 0, nframes*2-2 
    ang = ix*HO_twopi/float(nframes-1)
 
    if (ix <= nframes-1) then 
      call Healpix_GetRotation(R, ang, 0.d0, 0.d0)
    else
      call Healpix_GetRotation(R, ang-HO_twopi,HO_pi/2,0.d0)
    end if
!    call HealpixVis_map2ppm2(M,ppm,HealpixVis_proj_orth,1,400,400, R, mapix = pol)
     write (frame,"(I4.4)") ix

    call HealpixVis_ppm_write(ppm, trim(fileroot)//trim(frame)//'.ppm')    

  end do
  call HealpixVis_ppm_Free(ppm)

  end subroutine HealpixVis_map2Anim


 subroutine HealpixVis_Map2ppmfile(M, fname, projection, plot, n, symmetric, size_scale)
    type(HealpixMap), intent(in) :: M
    type(HealpixPPM):: ppm
    character(LEN=*), intent(in) :: fname
    character(LEN=*), intent(in), optional :: projection
    integer, intent(in), optional :: plot, n
    logical, intent(in), optional ::  symmetric
    real(SP), intent(in), optional :: size_scale
    integer splot, an
    logical symm
    real(sp) :: scale

    if (present(size_scale)) then
     scale= size_scale
    else
     scale=1
    end if 
    
    splot = splot_none
    if (present(plot)) then
      splot = plot
    end if

    if (present(n)) then
       an = n
    else
       an = 1
    end if
    
    if (present(symmetric)) then
     symm = symmetric
    else 
     symm = .false. 
    end if

    if (present(projection)) then 
     call HealpixVis_Map2PPM(M,ppm,projection,n=an,plot = splot, symmetric=symm, size_scale=scale)
    else 
     call HealpixVis_Map2PPM(M,ppm,n=an,plot = splot, symmetric=symm, size_scale=scale)
    end if

    call HealpixVis_ppm_write(ppm, fname)
    call HealpixVis_ppm_Free(ppm)

 end  subroutine HealpixVis_Map2ppmfile


 subroutine HealpixVis_MatlabSpinPlot(M, fname,n, nxx, nyy)
    use pix_tools, only : Ang2vec
    type(HealpixMap), intent(in) :: M
    character(len=*), intent(in) :: fname
    real ::  n2
    real(dp) :: theta, phi
    integer :: i, j, k, l, xsize, ysize, nxx, nyy
    integer, intent(in) :: n
    real :: x, y, xbase, ybase, xpix, ypix, xsubpix, ysubpix
    real(dp) vx,vy, ang, len
    integer :: pnum
    logical :: novec, inmap
    complex(dp) av


    open(1, file=fname, action='write', status='replace',form='formatted')
   ! write (1,3'(1I6)') M%spin
     
    xsize = nxx
    ysize = nyy
    n2 = real(n**2)
    xpix = 2.0/real(xsize)
    ypix = 2.0/real(ysize)
    xsubpix = xpix/real(2*n)
    ysubpix = ypix/real(2*n)
    do j = 0, ysize-1
      ybase = 1.0-real(j)*ypix
      do i = 0, xsize-1
        xbase = real(i)*xpix-1.0
        ! Inner loops sample map n^2 times in each pixel
        av=0
        novec = .false.
        do l = 0, n-1
          y = ybase - real(2*l+1)*ysubpix
          do k = 0, n-1
            x = xbase + real(2*k+1)*xsubpix
            call HealpixVis_proj_mol(x, y, inmap, theta, phi)
            if (.not.inmap) then
                   novec = .true.
            else
               pnum = HealpixMap_Ang2pix(M, theta, phi)
               
              if (real(M%SpinField(pnum)) == fmissval .or. M%SpinField(pnum) ==0) then
                 novec = .true.
              else
                 av = av + M%SpinField(pnum)
              end if
              
            end if
          end do
        end do
        !Matlab axes based at bottom left
        !we use ppm file or set(gca,'YDir','reverse') which sets positive being down and right
        !Q, U are e_theta and e_phi components
        !e_theta points down, and
        !we are looking at sky from inside, so e_phi points left, hence Q and -U
        if (.not. novec) then
         av = av /n2
          if (M%spin/=1) then
             vx = aimag(av)
             vy = real(av)          
             len = abs(av)
             ang = atan2(vx,vy)
             vx = -len*sin(ang/M%spin)
             vy = len*cos(ang/M%spin)
          else
           vx = -aimag(av)
           vy = real(av)
          end if 
          write (1,'(2I6,2e15.5)') i+1,j+1, vx, vy
        end if
       
      end do
    end do

   close(1)

 end subroutine HealpixVis_MatlabSpinPlot



  subroutine HealpixVis_map2ppm2(M, ppm, projection, n, nxx, nyy, pixarr, R, symmetric)
    use pix_tools, only : Ang2vec
    type(HealpixMap), intent(in) :: M
    real(sp), intent(in), target, optional :: pixarr(0:)
    integer, intent(in) :: n
    logical, intent(in), optional :: symmetric
    type(HealpixPPM), intent(inout) :: ppm
    real(dp), intent(in), optional :: R(3,3)
    !character(LEN=*), intent(in) :: colmap
    interface
      subroutine projection(x, y, inmap, theta, phi)
        integer, parameter :: dp = kind(1.0d0)
        real, intent(in) :: x, y
        logical, intent(out) :: inmap
        real(dp), intent(out) :: theta, phi
      end subroutine projection
    end interface

    real :: maprange, mapmin, mapmax, value, n2
    real(dp) :: theta, phi
    integer :: off,i, j, k, l, xsize, ysize, nxx, nyy
    real :: x, y, xbase, ybase, xpix, ypix, xsubpix, ysubpix
    real :: colour(3)
    real(sp), dimension(:), pointer :: PPix
    real(dp) vec(3)
    integer :: pnum
    logical :: inmap
    type(HealpixCMap) :: cmap

!!ppm starts at top left of image

    off = 0
    if (present(pixarr)) then
      PPix => pixarr
    else
      if (M%nmaps ==0) then
       PPix => M%Phi      
      else
       PPix => M%TQU(:,1)
       off = 1      
      end if
    end if

    call HealpixVis_cmap_rainbow(cmap)
    call HealpixVis_ppm_init(ppm, nxx, nyy)

    xsize = ppm%nx
    ysize = ppm%ny
    mapmin = minval(PPix)
    mapmax = maxval(PPix)
    if (present(symmetric)) then
     if (symmetric) then
      !symmetric about zero
       mapmin = min(mapmin,-mapmax)
       mapmax = -mapmin
     end if
    end if
   if (HealpixVis_force_range) then
     mapmin = -HealpixVis_abs_max
     maprange = 2*HealpixVis_abs_max
   else
     maprange = mapmax-mapmin
   end if
   if (VisDebugMsgs) then
     write (*,*) 'Map min = ',mapmin, 'Map max = ',maprange+mapmin
   end if

    n2 = real(n**2)
    xpix = 2.0/real(xsize)
    ypix = 2.0/real(ysize)
    xsubpix = xpix/real(2*n)
    ysubpix = ypix/real(2*n)
    do j = 0, ysize-1
      ybase = 1.0-real(j)*ypix
      do i = 0, xsize-1
        xbase = real(i)*xpix-1.0
        ! Inner loops sample map n^2 times in each pixel
        colour = 0
        do l = 0, n-1
          y = ybase - real(2*l+1)*ysubpix
          do k = 0, n-1
            x = xbase + real(2*k+1)*xsubpix
            call projection(x, y, inmap, theta, phi)
            if (.not.inmap) then
              colour = colour + ppm_background
            else
              if (present(R)) then
               call Ang2vec(theta, phi, vec)
               vec = matmul(R,vec)
               pnum = HealpixMap_Vec2pix(M, vec)
  
              else
               pnum = HealpixMap_Ang2pix(M, theta, phi)
              end if
              !write(*,*) i, j, x, y, phi, theta, pnum
              value=PPix(pnum+off) 
              
              if (value == fmissval) then  !.or. value ==0) then
                 colour = colour + ppm_null
              else
              colour = colour + HealpixVis_getcol(cmap,(value-mapmin)/maprange)
              end if
            end if
          end do
        end do
        ! Colour is average over n^2 colour samples
        ppm%rgb(:, i+1, j+1) = nint(255.0*colour/n2)
        ! Inner loops end here
      end do
    end do

   call HealpixVis_cmap_free(cmap)

  end subroutine HealpixVis_map2ppm2


    subroutine HealpixVis_MapMask2ppm(WM,WeightMap, fname, range) 
     Type(HealpixMap) :: WM, WeightMap
     character(LEN=*) :: fname 
     type(HealpixPPM):: ppm, ppmmask
     real, intent(in), optional :: range
     integer i,j
     logical :: old_force
     old_force = HealpixVis_force_range
     HealpixVis_force_range = .true.
     if (present(range)) then
      HealpixVis_abs_max = range
     else
      HealpixVis_abs_max = 500.
     end if
     call HealpixVis_Map2ppm(WM,ppm, n=4,plot = splot_none,symmetric = .true.)
     HealpixVis_force_range = .false.
     call HealpixVis_Map2ppm(WeightMap,ppmmask, n=1,plot = splot_none,symmetric = .false.)    
     do i=1,HealpixVis_DEF_WIDTH
      do j=1,HealpixVis_DEF_HEIGHT
       if (ppmmask%rgb(3,i,j)/=0 .and.  ppmmask%rgb(1,i,j)==0) then
         ppm%rgb(:,i,j)=0
       end if  
      end do
     end do 
    call HealpixVis_ppm_write(ppm,  fname)
    call HealpixVis_ppm_Free(ppm)
    call HealpixVis_ppm_Free(ppmmask)
    HealpixVis_force_range = old_force
    
    end subroutine HealpixVis_MapMask2ppm


  !======================================================================



  !==================================================================

  ! Equidistant cylindrical projection 

  subroutine HealpixVis_proj_ecp(x, y, inmap, theta, phi)

    ! Equidistant cylindrical projection

    real, intent(in) :: x, y
    logical, intent(out) :: inmap
    real(dp), intent(out) :: theta, phi

    stop 'need to check phi sign'
    inmap = .true.
    phi = x*HO_pi
    if (phi < 0.0_dp) phi = phi + HO_twopi 
    theta = (1.0_dp-y) * HO_pi/2.0_dp

  end subroutine HealpixVis_proj_ecp

  !==================================================================

  ! orthographic projection

  subroutine HealpixVis_proj_orth(x, y, inmap, theta, phi)
    use pix_tools, only: Vec2Ang
    real, intent(in) :: x, y
    logical, intent(out) :: inmap
    real(dp), intent(out) :: theta, phi
    real(dp) vec(3), r2    

    stop 'need to check phi sign'
 
    r2=  x**2 + y**2
    inmap = r2 <= 1
    if (inmap) then
       vec(1)=x
       vec(2)=y
       vec(3)=sqrt(1-r2)
       call Vec2Ang(vec,theta,phi)
    end if

  end subroutine HealpixVis_proj_orth


  !==================================================================

  ! Equal-area cylindrical projection

  subroutine HealpixVis_proj_cyl(x, y, inmap, theta, phi)

    real, intent(in) :: x, y
    logical, intent(out) :: inmap
    real(dp), intent(out) :: theta, phi

    stop 'need to check phi sign'

    inmap = .true.
    phi = x*HO_pi
    if (phi < 0.0_dp) phi = phi + HO_twopi
    theta = acos(real(y,dp))

  end subroutine HealpixVis_proj_cyl

  !==================================================================

  ! Sinusoidal projection

  subroutine HealpixVis_proj_sin(x, y, inmap, theta, phi)

    real, intent(in) :: x, y
    logical, intent(out) :: inmap
    real(dp), intent(out) :: theta, phi
    real(dp) :: snth

    stop 'need to check phi sign'

    theta = (1.0_dp-y) * HO_pi/2.0_dp
    snth = sin(theta)
    if (abs(x) > snth) then
      inmap = .false.
    else
      inmap = .true.
      phi = x*HO_pi/snth
      if (phi < 0.0) phi = phi + HO_twopi
    end if

  end subroutine HealpixVis_proj_sin

  !==================================================================

  ! Mollweide projection
  ! Looking out from inside so phi increases leftwards

  subroutine HealpixVis_proj_mol(x, y, inmap, theta, phi)

    real, intent(in) :: x, y
    logical,intent(out) :: inmap
    real(dp), intent(out) :: theta, phi
    real :: td
    real(dp) :: csth
    td = asin(y)
    csth = cos(td)
    if (abs(x) > csth) then
      inmap = .false.
    else
      inmap = .true.
      theta = HO_pi/2.0 - asin((2.0*td+sin(2.0*td))/HO_pi)
      phi = - x * HO_pi/csth
      if (phi < 0.0) phi = phi + HO_twopi
    end if

  end subroutine HealpixVis_proj_mol

  !==================================================================



  subroutine HealpixVis_cmap_init(cmap, n)
    type(HealpixCMap), intent(inout) :: cmap
    integer, intent(in) :: n
    integer err
    cmap%n = n
    deallocate(cmap%x,stat = err) 
    allocate(cmap%x(n), cmap%rgb(3,n))    
  end subroutine HealpixVis_cmap_init

  !======================================================================

  ! HealpixCMap destructor

  subroutine HealpixVis_cmap_free(cmap)
    type(HealpixCMap), intent(inout) :: cmap
    cmap%n = 0
    deallocate(cmap%x, cmap%rgb)
    nullify(cmap%x)
  end subroutine HealpixVis_cmap_free

  !======================================================================

  function HealpixVis_getcol(cmap, x) result(colour)

    type(HealpixCMap), intent(in) :: cmap
    real, intent(in) :: x
    real, dimension(3) :: colour
    integer :: i

    if (x < cmap%x(1)) then
      colour = cmap%rgb(:, 1)
    else if (x >= cmap%x(cmap%n)) then
      colour = cmap%rgb(:, cmap%n)
    else
      do i = 2, cmap%n
        if (x<cmap%x(i)) then
          colour = (cmap%rgb(:,i)*(x-cmap%x(i-1)) + &
              cmap%rgb(:,i-1)*(cmap%x(i)-x)) &
              / (cmap%x(i)-cmap%x(i-1))
          exit
        end if
      end do
    end if
  end function HealpixVis_getcol

  !======================================================================

  ! Allocate and set up greyscale colormap.

  subroutine HealpixVis_cmap_grey(cmap)
    type(HealpixCMap), intent(inout) :: cmap
    call HealpixVis_cmap_init(cmap, 2)
    cmap%x(1) = 0.0; cmap%rgb(:, 1) = (/0.0, 0.0, 0.0/)
    cmap%x(2) = 1.0; cmap%rgb(:, 2) = (/1.0, 1.0, 1.0/)
  end subroutine HealpixVis_cmap_grey

  !======================================================================

  ! Allocate and set up "angry orange" colourmap. 

  subroutine HealpixVis_cmap_orange(cmap)
    type(HealpixCMap), intent(inout) :: cmap
    call HealpixVis_cmap_init(cmap, 5)
    cmap%x(1) = 0.00; cmap%rgb(:, 1) = (/0.0, 0.0, 0.0/)
    cmap%x(2) = 0.25; cmap%rgb(:, 2) = (/0.5, 0.0, 0.0/)
    cmap%x(3) = 0.50; cmap%rgb(:, 3) = (/1.0, 0.5, 0.0/)
    cmap%x(4) = 0.75; cmap%rgb(:, 4) = (/1.0, 1.0, 0.5/)
    cmap%x(5) = 1.00; cmap%rgb(:, 5) = (/1.0, 1.0, 1.0/)

  end subroutine HealpixVis_cmap_orange

  !======================================================================

  ! Allocate and set up "rainbow" colourmap

  subroutine HealpixVis_cmap_rainbow(cmap)
    type(HealpixCMap), intent(inout) :: cmap
    call HealpixVis_cmap_init(cmap, 6)
    cmap%x(1) = 0.000; cmap%rgb(:, 1) = (/0.0, 0.0, 0.5/)
    cmap%x(2) = 0.125; cmap%rgb(:, 2) = (/0.0, 0.0, 1.0/)
    cmap%x(3) = 0.375; cmap%rgb(:, 3) = (/0.0, 1.0, 1.0/)
    cmap%x(4) = 0.625; cmap%rgb(:, 4) = (/1.0, 1.0, 0.0/)
    cmap%x(5) = 0.875; cmap%rgb(:, 5) = (/1.0, 0.0, 0.0/)
    cmap%x(6) = 1.000; cmap%rgb(:, 6) = (/0.5, 0.0, 0.0/)
  end subroutine HealpixVis_cmap_rainbow

  !======================================================================


end module HealpixVis
