subroutine gaussian_alm_1darray(nx,ny,D,lmin,lmax,alm,Cl,fix)
! Generate 1D-array random gaussian fields in 2D Fourier space for a given isotropic spectrum
  implicit none

  ![input]
  ! ny --- x and y grid number
  ! D(2)  --- x and y length
  ! iL(2) --- min/max multipoles of the random gaussian fields
  ! Cl(:) --- power spectrum
  integer, intent(in) :: nx, ny, lmin, lmax
  double precision, intent(in) :: Cl(0:lmax), D(2)

  ![in-output]
  ! alm(:) --- nx x ny size random gaussian fileds
  double complex, intent(out) :: alm(nx*ny)

  !(optional)
  ! fix   --- use sqrt{Cl} for alm (not random)
  logical, intent(in), optional :: fix

  !internal
  integer :: i, j, n, npixc, npix, l0, l1
  double precision :: x, y, dd(1:2), d0, l
  double precision, allocatable :: els(:), amp(:), amp2d(:)

  ! check
  if(size(iL)/=2)   stop 'error (gaussian_alm) : size of iL is not 2'
  if(mod(nx,2)/=0.or.mod(ny,2)/=0) stop 'error (gaussian_alm) : nx and/or ny should be even integers'
 
  call InitRandom(-1)

  npix = nx*ny

  !* make cl on 2d grid
  allocate(amp(size(Cl)),els(npix),amp2d(npix))
  amp = 0d0
  d0 = D(1)*D(2)
  do i = iL(1), iL(2)
    if (Cl(i)<0d0) stop 'error: cl is negative'
    amp(i) = dsqrt(d0*Cl(i)*0.5d0)  ! \sqrt(\delta(l=0)*Cl(l)/2)
  end do
  if(present(fix)) amp = dsqrt(d0*Cl)
  els = elarray(nn,D)
  amp2d = 0d0
  do n = 1, npix
    if(iL(1)>els(n).or.els(n)>iL(2)-1) cycle
    l0 = int(els(n))
    l1 = l0 + 1
    amp2d(n) = amp(l0) + (els(n)-l0)*(amp(l1)-amp(l0))
  end do
  deallocate(amp,els)

  !* alm=0 if i=1 or j=1 for symmetry
  !* center: (ic,jc) = (nx/2+1, ny/2+1)
  alm = 0d0
  npixc = ny*nx/2 + ny/2 + 1

  !* dx, dy in l-space
  dd = twopi/D

  !* check
  if(dd(1)*(nx*0.5-1)<iL(2).or.dd(2)*(ny*0.5-1)<iL(2)) then
    write(*,*) 'error: inclusion of Fourier mode is incorrect'
    write(*,*) 'maximum ell should be lower than', dd(1)*(nx*0.5-1), 'or', dd(2)*(ny*0.5-1)
    stop
  end if

  ! half side (i < ic)
  n = ny+1 ! adding from j=1 to ny
  do i = 2, nx/2
    x = dd(1)*dble(i-1-nx*0.5d0)
    n = n + 1 ! j=1
    do j = 2, ny
      y = dd(2)*dble(j-1-ny*0.5d0)
      l = dsqrt(x**2+y**2)
      if(l>=iL(1).and.l<=iL(2)) then
        if(present(fix)) then
          alm(n) = amp2d(n)
        else
          alm(n) = cmplx(Gaussian1(),Gaussian1())*amp2d(n)
        end if
      end if
      n = n + 1
    end do
  end do

  ! values on axis (i=ic) but avoid the origin (ell=0) 
  i = nx/2+1
  ! x=0
  n = n + 1  ! j = 1
  do j = 2, ny/2
    y = dd(2)*dble(j-1-ny*0.5d0)
    l = abs(y)
    if(l>=iL(1).and.l<=iL(2)) then
      if(present(fix)) then
        alm(n) = amp2d(n)
      else
        alm(n) = cmplx(Gaussian1(),Gaussian1())*amp2d(n)
      end if
    end if
    n = n + 1
  end do

  !complex conjugate
  do n = npixc+1, nx*ny
    alm(n) = conjg(alm(2*npixc-n))
  end do

  deallocate(amp2d)
  

end subroutine gaussian_alm_1darray


subroutine maskpoint_generate(nn,D,Nm,px,py)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), Nm
  double precision, intent(in) :: D(2)
  double precision, intent(out) :: px(:), py(:)
  !internal
  integer :: i, j, n, m

  m = Nm
  write(*,*) 'generate masked points'
  call InitRandom(-1)
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      if(m>0.or.dble(m)/dble(nn(1)*nn(2)-n)>ranmar()) then
        px(m) = dble(i-1-(nn(1)-1)*0.5)*D(1)/dble(nn(1))
        py(m) = dble(j-1-(nn(2)-1)*0.5)*D(2)/dble(nn(2))
        m = m - 1
      end if
      n = n + 1
    end do
  end do

end subroutine maskpoint_generate


function Mwin(mr,t,mfac) 
!* apodized mask window function
  implicit none
  !I/O
  double precision, intent(in) :: mr, t
  double precision, intent(in), optional :: mfac
  !internal
  double precision :: t0, Mwin, x, b

  b = mr
  if(present(mfac)) b = b*mfac

  if (t<mr) then
    Mwin = 1d0
  else if (t>=mr.and.t<b.and..not.b==mr) then
    x = (b-t)/(b-mr)
    Mwin = x - sin(2*pi*x)/(2*pi)
  else 
    Mwin = 0d0
  end if

end function Mwin


subroutine mask_generate(nn,D,px,py,mr,M)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), mr, px(:), py(:)
  double precision, intent(inout) :: M(:)
  !internal
  integer :: i, j, p, n
  double precision :: xi, xj, dx, dy, sx ,sy

  sx = D(1)/dble(nn(1))
  sy = D(2)/dble(nn(2))

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      do p = 1, size(px)
        xi = dble(i-1-(nn(1)-1)*0.5)*sx
        xj = dble(j-1-(nn(2)-1)*0.5)*sy
        dx = abs(xi-px(p))
        dy = abs(xj-py(p))
        M(n) = (1d0-Mwin(mr*sx,dx)*Mwin(mr*sy,dy))*M(n)
      end do
      n = n + 1
    end do
  end do

end subroutine mask_generate


subroutine window_generate(nn,D,W,mr,mn,f,ap,cut)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: nn(2)
  integer, intent(in), optional :: mn
  double precision, intent(in) :: D(2)
  double precision, intent(in), optional :: mr, ap, cut
  double precision, intent(inout) :: W(:)
  !internal
  integer :: i, j, n, maskn
  double precision :: a, c, xi, xj, dx, dy, sx ,sy
  double precision, allocatable :: px(:), py(:), M(:)

  sx = D(1)/dble(nn(1))
  sy = D(2)/dble(nn(2))

  a  = 1d0
  if (present(ap))  a = ap

  c  = 1d0
  if (present(cut)) c = cut

  call window_sin(nn,D,W,a,c)

  allocate(M(nn(1)*nn(2)));  M = 1d0
  if(present(mr).and..not.mr==0) then
    if(present(f)) then
      maskn = filelines(f)
      allocate(px(maskn),py(maskn))
      write(*,*) 'read masked points from a file, mask points =', maskn
      call loadtxt(f,px,py)
    else if(present(mn)) then
      maskn = mn
      allocate(px(maskn),py(maskn))
      call maskpoint_generate(nn,D,maskn,px,py)
      call savetxt('genpoints.dat',px,py)
    else
      stop 'error: require number of masks to be generated'
    end if
    call mask_generate(nn,D,px,py,mr,M)
    deallocate(px,py)
  end if

  W = W*M

  deallocate(M)

end subroutine window_generate

subroutine gaussbeam_2d(t,els,beam,eL)
!* gaussain beam in 2D
  implicit none
  double precision, intent(in) :: els(:), t
  double precision, intent(out) :: beam(:)
  integer, intent(in), optional :: eL(1:2)
  !internal
  integer :: n, cL(2)

  cL = [0,int(maxval(abs(els)))+1]
  if (present(eL)) cL = eL

  do n = 1, size(els)
    if (cL(1)<=els(n).and.els(n)<=cL(2)) beam(n) = dexp(-els(n)**2*t**2/16d0/dlog(2d0))
  end do

end subroutine gaussbeam_2d


subroutine fourier_filter_1d(els,eL,alm,Fl,Fn)
! * assign zerros outside of annular region
! * optionally multiply a filter function
  implicit none
  integer, intent(in) :: eL(1:2)
  double precision, intent(in) :: els(:)
  double complex, intent(inout) :: alm(:)
  double precision, intent(in), optional :: Fl(:)
  double precision, intent(in), optional :: Fn(:)
! [internal]
  integer :: n, l

  if (present(Fl).and.present(Fn)) stop 'error: two filters appeared (fourier filtering)'

  do n = 1, size(alm)
    l = int(els(n))
    if (l>=eL(1).and.l<=eL(2)) then
      if (present(Fl)) alm(n) = alm(n)*Fl(l)
      if (present(Fn)) alm(n) = alm(n)*Fn(n)
    else
      alm(n) = 0d0
    end if
  end do

end subroutine fourier_filter_1d


subroutine fourier_filter_2d(els,eL,alm,Fl,Fn)
! * same as 'fourier_filter_1d' but for 2d array
  implicit none
  integer, intent(in) :: eL(1:2)
  double precision, intent(in) :: els(:)
  double complex, intent(inout) :: alm(:,:)
  double precision, intent(in), optional :: Fl(:,:)
  double precision, intent(in), optional :: Fn(:,:)
! [internal]
  integer :: n, l

  if (present(Fl).and.present(Fn)) stop 'error: two filters appeared (fourier filtering)'

  do n = 1, size(alm,dim=2)
    l = int(els(n))
    if (l>=eL(1).and.l<=eL(2)) then
      if (present(Fl)) alm(:,n) = alm(:,n)*Fl(:,l)
      if (present(Fn)) alm(:,n) = alm(:,n)*Fn(:,n)
    else
      alm(:,n) = 0d0
    end if
  end do

end subroutine fourier_filter_2d


subroutine elarray_wcut(nn,D,oL,els,lfac,elx,linv)
  implicit none
  integer, intent(in) :: nn(1:2), oL(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(out), optional :: els(:), lfac(:), elx(:), linv(:)
  integer :: n, npix
  double precision, allocatable :: els0(:), lfac0(:), linv0(:)

  npix = nn(1)*nn(2)
  allocate(els0(npix),lfac0(npix),linv0(npix)); lfac0=0d0; linv0=0d0
  els0 = elarray(nn,D)
  do n = 1, npix
    if (oL(1)<=els0(n).and.els0(n)<=oL(2)) lfac0(n) = 2d0/els0(n)**2
    if (oL(1)<=els0(n).and.els0(n)<=oL(2)) linv0(n) = 1d0/els0(n)
  end do
  if (present(els))  els  = els0
  if (present(lfac)) lfac = lfac0
  if (present(elx))  elx  = elarray_x(nn,D)
  if (present(linv)) linv = linv0
  deallocate(els0,lfac0,linv0)

end subroutine elarray_wcut


subroutine array_fine(nn,s,arr1,arr2,f) 
!* make finer array with a same value
  implicit none
!
! [input]
! nn   --- x and y grids
! s    --- increasing factors 
! arr1 --- finer array with nn(1) x nn(2) grids
  integer, intent(in) :: nn(1:2), s(1:2)
  double precision, intent(in) :: arr1(:)
!
! [output]
! arr2 --- finer array with nn(1)*s(1) x nn(2)*s(2) grids
  double precision, intent(out) :: arr2(:)
!
! (optional)
! f    --- file names for input and output arrays
  character(*), intent(in), optional :: f(1:2)
!
! [internal]
  integer :: mm(2), i, j
  double precision, allocatable :: map1(:,:), map2(:,:)

  mm = nn*s
  allocate(map1(nn(1),nn(2)),map2(mm(1),mm(2)))
  if (present(f))  call savetxt(f(1),arr1)

  map1 = reshape(arr1,nn,order=[2,1])
  do i = 1, nn(1)
    do j = 1, nn(2)
      map2((i-1)*s(1)+1:(i-1)*s(1)+s(1),(j-1)*s(2)+1:(j-1)*s(2)+s(2)) = map1(i,j)
    end do
  end do
  arr2 = reshape(transpose(map2),[mm(1)*mm(2)])
  if (present(f))  call savetxt(f(2),arr2)

  deallocate(map1,map2)

end subroutine array_fine


subroutine map_smoothing_2D(nn,rn,s,map,smap)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double precision, intent(in) :: map(:,:)
  double precision, intent(out) :: smap(:,:)
  !internal
  integer :: i, j

  do i = 1, rn(1)
    do j = 1, rn(2)
      smap(i,j) = sum(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))/dble(s**2)
    end do
  end do

end subroutine map_smoothing_2D


subroutine map_smoothing_2D_adv(nn,rn,s,map,smap)
!* An extention of map_smoothing_2D
!* set zero if one of pixels has zero value
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double precision, intent(in) :: map(:,:)
  double precision, intent(out) :: smap(:,:)
  !internal
  integer :: i, j

  do i = 1, rn(1)
    do j = 1, rn(2)
      if(product(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))==0) then
        smap(i,j) = 0d0
      else
        smap(i,j) = sum(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))/dble(s**2)
      end if
    end do
  end do

end subroutine map_smoothing_2D_adv


subroutine map_smoothing_1D(nn,rn,s,map1D,smap1D)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double precision, intent(in) :: map1D(:)
  double precision, intent(out) :: smap1D(:)
  !internal
  double precision :: map(nn(1),nn(2)), smap(rn(1),rn(2))

  map = reshape(map1D,[nn(1),nn(2)],order=[2,1])
  call map_smoothing_2D(nn,rn,s,map,smap)
  smap1D = reshape(transpose(smap),[rn(1)*rn(2)])

end subroutine map_smoothing_1D


subroutine map_smoothing_1D_cmplx(nn,rn,s,map1D,smap1D)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double complex, intent(in) :: map1D(:)
  double complex, intent(out) :: smap1D(:)
  !internal
  double complex :: map(nn(1),nn(2)), smap(rn(1),rn(2))

  map = reshape(map1D,[nn(1),nn(2)],order=[2,1])
  call map_smoothing_2D_cmplx(nn,rn,s,map,smap)
  smap1D = reshape(transpose(smap),[rn(1)*rn(2)])

end subroutine map_smoothing_1D_cmplx


subroutine map_smoothing_2D_cmplx(nn,rn,s,map,smap)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double complex, intent(in) :: map(:,:)
  double complex, intent(out) :: smap(:,:)
  !internal
  integer :: i, j

  do i = 1, rn(1)
    do j = 1, rn(2)
      smap(i,j) = sum(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))/dble(s**2)
    end do
  end do

end subroutine map_smoothing_2D_cmplx



