!////////////////////////////////////////////////////!
! * Utils for 2D Grid Analysis
!////////////////////////////////////////////////////!

module utils
  use random,    only: InitRandom, ranmar, Gaussian1
  use constants, only: iu, pi
  use general,   only: check_positive
  use grid2d,    only: elxy, elarray, elarrays_1d, elarrays_2d, wap, make_lmask
  use pstool,    only: binned_ells, power_binning, cb2cl
  implicit none

  private check_positive
  private InitRandom, ranmar, gaussian1
  private iu, pi
  private elxy, elarray, elarrays_1d, elarrays_2d, wap, make_lmask
  private binned_ells, power_binning, cb2cl


contains 

!//// Fourier modes ////!

subroutine el2d(nx,ny,D,els)
!*  Return absolute value of multipole in 2D grids
!* 
!*  Args:
!*    :nx, ny (int):  number of Lx and Ly grids
!*    :D[xy] (double):  map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*  
!*  Returns:
!*    :els[nx,ny] (double):  absolute value of Fourier mode, (Lx**2+Ly**2)**0.5, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(out), dimension(nx,ny) :: els
  !internal
  integer :: i, j
  double precision :: lx, ly

  do i = 1, nx
    lx = elxy(i,nx,D(1))
    do j = 1, ny
      ly = elxy(j,ny,D(2))
      if (lx==0.and.ly==0) cycle
      els(i,j) = dsqrt(lx**2+ly**2)
    end do
  end do

end subroutine el2d


subroutine elarrays(nx,ny,D,elx,ely,els,eli)
!*  Return Lx, Ly, absolute value of multipole, and its inverse in 2D grids
!* 
!*  Args:
!*    :nx, ny (int):  number of Lx and Ly grids
!*    :D[xy] (double):  map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*  
!*  Returns:
!*    :elx[nx,ny] (double) : Lx, with bounds (nx,ny)
!*    :ely[nx,ny] (double) : Ly, with bounds (nx,ny)
!*    :els[nx,ny] (double) : absolute value of Fourier mode, (Lx**2+Ly**2)**0.5, with bounds (nx,ny)
!*    :eli[nx,ny] (double) : inverse of els, with bounds (nx,ny)
!*
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(out), dimension(nx,ny) :: elx, ely, els, eli

  call elarrays_2d((/nx,ny/),D,elx,ely,els,eli)

end subroutine elarrays


subroutine elmask(nx,ny,D,lmask,lmin,lmax,lxcut,lycut)
!* Return mask in 2D Fourier space. The lmask is unity at lmin<=|L|<=lmax, |Lx|>=lxcut, |Ly|>=lycut, and otherwize zero. 
!*
!* Args: 
!*    :nx, ny (int):  number of Lx and Ly grids
!*    :D[xy] (double):  map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*
!* Args(optional):
!*    :lmin/lmax (int)  : Minimum/Maximum of multipoles
!*    :lxcut/lycut (int): Remove |Lx|<lxcut / |Ly|<lycut cut of multipoles
!*  
!*  Returns:
!*    :lmask[nx,ny] (double) : Mask, with bounds (nx,ny)
!*
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(out), dimension(nx,ny) :: lmask
  double precision, intent(in), optional :: lmin, lmax, lxcut, lycut
  !f2py integer :: lmin  = 0
  !f2py integer :: lmax  = 1000
  !f2py integer :: lxcut = 0
  !f2py integer :: lycut = 0
  integer :: i, j
  double precision :: els(nx,ny), elx(nx,ny), ely(nx,ny)

  call elarrays_2d((/nx,ny/),D,elx,ely,els)

  lmask = 1d0
  do i = 1, nx
    do j = 1, ny
      if (present(lmin).and.els(i,j)<lmin)        lmask(i,j) = 0d0
      if (present(lmax).and.els(i,j)>lmax)        lmask(i,j) = 0d0
      if (present(lxcut).and.abs(elx(i,j))<lxcut) lmask(i,j) = 0d0
      if (present(lycut).and.abs(ely(i,j))<lycut) lmask(i,j) = 0d0
    end do
  end do

end subroutine elmask

!//// Compute Cl from Fourier modes ////!

subroutine alm2bcl(bn,oL,nx,ny,D,Cb,alm1,alm2,spc)
!*  Compute angular power spectrum from Fourier modes, with multipole binning
!* 
!*  Args:
!*    :bn (int)      : number of multipole bin
!*    :oL[2] (int)   : minimum and maximum multipoles of the output cl
!*    :nx, ny (int)  : number of Lx and Ly grids
!*    :D[xy] (double): map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :alm1[nx,ny] (dcmplx): Fourier mode, with bounds (nx,ny)
!* 
!*  Args(optional):
!*    :alm2[nx,ny] (dcmplx): Fourier mode, with bounds (nx,ny), default to None
!*    :spc (str)           : type of multipole binning, i.e., linear spacing (spc='', default), or log spacing (spc='log')
!*
!*  Returns:
!*    :Cb[bin] (double) : angular power spectrum with multipole binning, with bounds (bn)
!*
  implicit none
  !inputs
  integer, intent(in) :: bn, nx, ny
  integer, intent(in), dimension(2) :: oL
  double precision, intent(in), dimension(2) :: D
  double complex, intent(in), dimension(nx,ny) :: alm1
  !optional
  character(*), intent(in), optional :: spc
  double complex, intent(in), dimension(nx,ny), optional :: alm2
  !f2py character(*) :: spc=''
  !f2py double complex :: alm2 = 0
  !docstr :: alm2 = alm1
  !outputs
  double precision, intent(out), dimension(bn) :: Cb
  !internal
  character(8) :: spc0
  double precision, allocatable :: C(:,:)
  double complex, allocatable :: alm0(:,:)

  spc0  = ''
  if (present(spc))  spc0 = spc

  !* 2D power spectrum
  allocate(C(nx,ny),alm0(nx,ny))
  alm0 = alm1
  if (present(alm2).and.sum(abs(alm2))==0) alm0 = alm1
  if (present(alm2).and.sum(abs(alm2))/=0) alm0 = alm2
  C = (alm1*conjg(alm0)+alm0*conjg(alm1))*0.5d0/(D(1)*D(2)) 
  deallocate(alm0)

  !* to 1D power spectrum
  call c2d2bcl(nx,ny,D,C,bn,oL,Cb,spc=spc0)

  deallocate(C)


end subroutine alm2bcl


subroutine c2d2bcl(nx,ny,D,c2d,bn,oL,Cb,spc)
!*  Return 1D angular power spectrum with multipole binning from a 2D power spectrum
!*
!*  Args:
!*    :nx, ny (int)       : number of Lx and Ly grids
!*    :D[xy] (double)     : map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :c2d[nx,ny] (double): 2D power spectrum, with bounds (nx,ny)
!*    :bn (int)           : number of multipole bin
!*    :oL[2] (int)        : minimum and maximum multipoles of the output cl
!*    
!*  Args(optional):
!*    :spc (str) : type of multipole binning, i.e., linear spacing (spc='', default), or log spacing (spc='log')
!*
!*  Returns:
!*    :Cb[bin] (double) : angular power spectrum with multipole binning, with bounds (bn)
!*
  implicit none
  !I/O
  integer, intent(in) :: bn, nx, ny
  integer, intent(in), dimension(2) :: oL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: c2d
  double precision, intent(out), dimension(bn) :: Cb
  !optional
  character(*), intent(in), optional :: spc
  !f2py character(*) :: spc=''
  !internal
  character(8) :: spc0
  double precision :: bp(bn+1), b(bn), els(nx,ny), lmask(nx,ny)
  double precision, dimension(bn) :: vAb

  spc0 = ''
  if (present(spc)) spc0 = spc

  call binned_ells(oL,bp,b,spc0)

  call make_lmask((/nx,ny/),D,oL,lmask)
  call el2d(nx,ny,D,els)
  call power_binning(bp,els,lmask*c2d,lmask,Cb,vAb)

end subroutine c2d2bcl


! Power spectrum interpolation

subroutine cl2c2d(nx,ny,D,lmin,lmax,Cl,c2d)
!*  Assign values of 1D angular power spectrum on to 2D grid with linear interpolation
!*
!*  Args: 
!*    :nx, ny (int)   : number of Lx and Ly grids
!*    :D[xy] (double) : map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :lmin (int)     : minimum multipole of cl to be interpolated
!*    :lmax (int)     : maximum multipole of cl to be interpolated
!*    :Cl[l] (double) : 1D power spectrum, with bounds (0:lmax)
!*
!*  Returns:
!*    :c2d[nx,ny] (double): 2D power spectrum, with bounds (nx,ny)
!* 
  implicit none
  !I/O
  integer, intent(in) :: nx, ny, lmin, lmax
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(0:lmax) :: Cl
  double precision, intent(out), dimension(nx,ny) :: c2d
  !internal
  logical :: p
  integer :: i, j, l0, l1
  double precision :: els(nx,ny)

  !check positivity
  call check_positive(Cl,p)

  call el2d(nx,ny,D,els)

  c2d = 0d0
  do i = 1, nx
    do j = 1, ny
      if (lmin>els(i,j).or.els(i,j)>lmax-1) cycle
      l0 = int(els(i,j))
      l1 = l0 + 1
      c2d(i,j) = Cl(l0) + (els(i,j)-l0)*(Cl(l1)-Cl(l0))
      if (c2d(i,j)>=0d0.or..not.p) cycle
      write(*,*) Cl(l0), Cl(l1), l0, els(i,j)
      stop 'error (cl2c2d): interpolated Cl is negative'
    end do
  end do

end subroutine cl2c2d


subroutine cb2c2d(bn,bc,nx,ny,D,lmin,lmax,Cb,C2d,method)
!*  Assign values of 1D angular power spectrum on to 2D grid with linear interpolation
!*
!*  Args: 
!*    :bn (int)        : number of multipole bins
!*    :bc[bin] (double): multipole bin center, with bounds (bn)
!*    :nx, ny (int)    : number of Lx and Ly grids
!*    :D[xy] (double)  : map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :lmin (int)      : minimum multipole of cl to be interpolated
!*    :lmax (int)      : maximum multipole of cl to be interpolated
!*    :Cb[bin] (double): 1D power spectrum with multipole binning, with bounds (bn)
!*
!*  Args(optional):
!*    :method (str) : interpolation method from binned to unbinned angular spectrum, i.e., spline ('', default), or linear ('linear') interpolation
!*
!*  Returns:
!*    :c2d[nx,ny] (double): 2D power spectrum, with bounds (nx,ny)
!* 
  implicit none
  !I/O
  integer, intent(in) :: bn, nx, ny, lmin, lmax
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(bn) :: bc, Cb
  double precision, intent(out), dimension(nx,ny) :: C2d
  !optional
  character(*), intent(in), optional :: method
  !f2py character(*) :: method=''
  !internal
  character(8) :: m
  double precision, allocatable :: Cl(:)

  allocate(Cl(0:lmax)); Cl=0d0

  m = ''
  if(present(method)) m   = method

  !interpolate Cb -> Cl
  call cb2cl(bc,Cb,Cl(1:lmax),method=m)

  !interpolate Cl -> C2d
  call cl2c2d(nx,ny,D,lmin,lmax,Cl,c2d)

  deallocate(Cl)

end subroutine cb2c2d


!//// Gaussian field generation ////!

subroutine gauss1alm(nx,ny,D,lmin,lmax,Cl,alm)
!*  Generate random gaussian fields in 2D Fourier space for a given isotropic spectrum, satisfying a^*_l = a_{-l}
!*
!*  Args:
!*    :nx, ny (int)    : number of Lx and Ly grids
!*    :D[xy] (double)  : map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :lmin (int)      : minimum multipole of cl to be interpolated
!*    :lmax (int)      : maximum multipole of cl to be interpolated
!*    :Cl[l] (double) : 1D power spectrum, with bounds (0:lmax)
!*
!*  Returns:
!*    :alm[lx,ly] (dcmplx): random gaussian fields on 2D Fourier plane, with bounds (nx,ny)
!*
  implicit none
  integer, intent(in) :: lmin, lmax, nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(0:lmax) :: Cl
  double complex, intent(out), dimension(nx,ny) :: alm
  !internal
  integer :: i, j, n, l0, l1
  integer :: ijmin(2), ijmax(2)
  double precision :: x, y, dx, dy, d0, l
  double precision, allocatable :: amp(:), amp2d(:,:)

  ! check
  if(mod(nx,2)/=0.or.mod(ny,2)/=0) write(*,*) 'WARNING (gauss1alm) : current code assumes nx, ny to be even'
 
  call initrandom(-1)

  ! make cl on 2d grid
  allocate(amp(size(Cl)),amp2d(nx,ny))
  amp = 0d0
  d0 = D(1)*D(2)
  do i = lmin, lmax
    if (Cl(i)<0d0) then
      write(*,*) 'error (gauss1alm): cl is negative', cl(i), i
      stop
    end if
    amp(i) = dsqrt(d0*Cl(i)*0.5d0)  ! \sqrt(\delta(l=0)*Cl(l)/2)
  end do
  amp2d = 0d0
  do i = 1, nx
    x = elxy(i,nx,D(1))
    do j = 1, ny
      y = elxy(j,ny,D(2))
      l = dsqrt(x**2+y**2)
      if(lmin>l.or.l>lmax-1) cycle
      l0 = int(l)
      l1 = l0 + 1
      amp2d(i,j) = amp(l0) + (l-l0)*(amp(l1)-amp(l0))
    end do
  end do
  deallocate(amp)

  ! alm=0 if i=1 or j=1 for symmetry
  ! center: (ic,jc) = (nx/2+1, ny/2+1)
  alm = 0d0

  ! maximum nn
  !ijmin = max ( 2, int( nn(:)*0.5d0 + 1 - iL(2)*D(:)/twopi ) )
  !ijmax = min ( nn(:), int( nn(:)*0.5d0 + 1 + iL(2)*D(:)/twopi ) + 1 )
  ijmin = 2
  ijmax = (/nx,ny/)

  ! check
  if(elxy(nx,nx,D(1))<lmax.or.elxy(ny,ny,D(2))<lmax) then
    write(*,*) 'error (gauss1alm): inclusion of Fourier mode is incorrect'
    write(*,*) 'maximum ell should be lower than', elxy(nx,nx,D(1)), 'or', elxy(ny,ny,D(2))
    stop
  end if

  ! half side (i < ic)
  do i = ijmin(1), nx/2
    x = elxy(i,nx,D(1))
    do j = ijmin(2), ijmax(2)
      y = elxy(j,ny,D(2))
      l = dsqrt(x**2+y**2)
      if(l<lmin.or.l>lmax) cycle
      alm(i,j) = cmplx(Gaussian1(),Gaussian1())*amp2d(i,j)
    end do
  end do

  ! values on axis (i=ic) but avoid the origin (ell=0) 
  i = nx/2+1
  ! x=0
  do j = ijmin(2), ny/2
    y = elxy(j,ny,D(2))
    l = abs(y)
    if(l<lmin.or.l>lmax) cycle
    alm(i,j) = cmplx(Gaussian1(),Gaussian1())*amp2d(i,j)
  end do
  do j = ny/2+2, ijmax(2)
    alm(i,j) = conjg(alm(i,ny-j+2))
  end do
  
  ! the other half side (i>ic)
  do i = nx/2+2, ijmax(1)
    do j = ijmin(2), ijmax(2)
      alm(i,j) = conjg(alm(nx-i+2,ny-j+2))
    end do
  end do

  deallocate(amp2d)


end subroutine gauss1alm


subroutine gauss2alm(nx,ny,D,lmin,lmax,TT,TE,EE,tlm,elm)
!*  Generate two correlated random gaussian fields in 2D Fourier space for a given isotropic spectrum
!*
!*  Args:
!*    :nx, ny (int)    : number of Lx and Ly grids
!*    :D[xy] (double)  : map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :lmin (int)      : minimum multipole of cl to be interpolated
!*    :lmax (int)      : maximum multipole of cl to be interpolated
!*    :TT[l] (double)  : the 1st 1D power spectrum, with bounds (0:lmax)
!*    :TE[l] (double)  : the cross 1D power spectrum, with bounds (0:lmax)
!*    :EE[l] (double)  : the 2nd 1D power spectrum, with bounds (0:lmax)
!*
!*  Returns:
!*    :tlm[lx,ly] (dcmplx): the 1st random gaussian fields on 2D Fourier plane, with bounds (nx,ny)
!*    :elm[lx,ly] (dcmplx): the 2nd random gaussian fields on 2D Fourier plane, with bounds (nx,ny)
!*
  implicit none
  integer, intent(in) :: nx, ny, lmin, lmax
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(0:lmax) :: TT, TE, EE
  double complex, intent(out), dimension(nx,ny) :: tlm, elm
  integer :: i, j
  double precision, allocatable :: uni(:), TT2d(:,:), TE2d(:,:), EE2d(:,:)
  double complex, allocatable :: ulm(:,:)

  allocate(uni(0:lmax),ulm(nx,ny)); uni=1d0; ulm=0d0

  call gauss1alm(nx,ny,D,lmin,lmax,TT,tlm)
  call gauss1alm(nx,ny,D,lmin,lmax,uni,ulm)

  allocate(TT2d(nx,ny),TE2d(nx,ny),EE2d(nx,ny))
  call cl2c2d(nx,ny,D,lmin,lmax,TT,TT2d)
  call cl2c2d(nx,ny,D,lmin,lmax,TE,TE2d)
  call cl2c2d(nx,ny,D,lmin,lmax,EE,EE2d)

  elm = 0d0
  do i = 1, nx
    do j = 1, ny
      if (TT2d(i,j)==0d0) cycle
      elm(i,j) = (TE2d(i,j)/TT2d(i,j))*tlm(i,j) + ulm(i,j)*dsqrt(EE2d(i,j)-TE2d(i,j)**2/TT2d(i,j))
    end do
  end do

  deallocate(TT2d,TE2d,EE2d,ulm,uni)

end subroutine gauss2alm


!//// Window function //////////////////////////////////////////////////////////////////////////////!

subroutine window_sin(nx,ny,D,W,ap,cut)
!*  Return a sin window function.
!*
!*  Args:
!*    :nx, ny (int)  : Number of Lx and Ly grids
!*    :D[xy] (double): Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*
!*  Args(Optional):
!*    :ap (double) : Apodization parameter defined by apodized-range = (1-ap) x (cut)mapsize, from 0 (full apodization) to 1 (no apodization). Default to 1.
!*    :cut (double): Map cut scale defined by cutmapsize = cut x mapsize, from 0 (full cut) to 1 (no cut). Default to 1.
!*
!*  Return:
!*    :W[x,y] (double): Window function, with bounds (nx,ny)
!*
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(out), dimension(nx,ny) :: W
  !optional
  double precision, intent(in), optional :: ap, cut
  !f2py double precision :: ap = 1
  !f2py double precision :: cut = 1
  !internal
  integer :: i, j
  double precision :: a, c, xi, xj, sx ,sy

  sx = D(1)/dble(nx)
  sy = D(2)/dble(ny)

  a  = 1d0
  if (present(ap)) a = ap

  c  = 1d0
  if (present(cut)) c = cut

  do i = 1, nx
    do j = 1, ny
      xi = dble(i-1-(nx-1)*0.5)*sx
      xj = dble(j-1-(ny-1)*0.5)*sy
      W(i,j) = wap(D(1),abs(xi),a,c)*wap(D(2),abs(xj),a,c)
    end do
  end do

end subroutine window_sin


subroutine window_norm(nx,ny,wind,num,wn)
  implicit none
  integer, intent(in) :: nx, ny, num
  double precision, intent(in), dimension(nx,ny) :: wind
  double precision, intent(out), dimension(0:num) :: wn
  integer :: n

  wn(0) = 1d0
  do n = 1, num
    wn(n) = sum(wind**n)/dble(nx*ny)
  end do

end subroutine window_norm


subroutine window_norm_x(nx,ny,W1,W2,num,Wn)
  implicit none
  integer, intent(in) :: nx, ny, num
  double precision, intent(in), dimension(nx,ny) :: W1, W2
  double precision, intent(out), dimension(0:num,0:num) :: Wn
  integer :: n, m

  do n = 0, num
    do m = 0, num
      Wn(n,m) = sum(W1**n*W2**m)/dble(nx*ny)
    end do
  end do

end subroutine window_norm_x


subroutine rotation(nx,ny,rot,QU,rQU,rtype)
  implicit none
  character(*), intent(in) :: rtype
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(nx,ny) :: rot
  double precision, intent(in), dimension(nx,ny,2) :: QU
  double precision, intent(out), dimension(nx,ny,2) :: rQU
  double precision, allocatable :: tmp(:,:,:)  

  allocate(tmp(nx,ny,2))
  select case(rtype)
  case('l')
    tmp(:,:,1) = QU(:,:,1) - QU(:,:,2)*2d0*rot
    tmp(:,:,2) = QU(:,:,2) + QU(:,:,1)*2d0*rot
  case('f')
    tmp(:,:,1) = QU(:,:,1)*dcos(2d0*dble(rot)) - QU(:,:,2)*dsin(2d0*dble(rot))
    tmp(:,:,2) = QU(:,:,2)*dcos(2d0*dble(rot)) + QU(:,:,1)*dsin(2d0*dble(rot))
  case default
    stop 'error (rotation): rotation type not specified'
  end select
  rQU = tmp
  deallocate(tmp)

end subroutine rotation


subroutine get_angle(nx,ny,D,theta,phi)
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in), dimension(2) :: D
  double precision, intent(out), dimension(ny) :: theta
  double precision, intent(out), dimension(nx) :: phi
  integer :: a, b

  do a = 1, nx
    phi(a)   = dble(a-1-(nx-1)*0.5)*D(1)/dble(nx)
  end do
  do b = 1, ny
    theta(b) = pi/2d0 - dble(b-1-(ny-1)*0.5)*D(2)/dble(ny)
  end do

end subroutine get_angle


subroutine cutmap(ox,oy,cx,cy,omap,cmap)
  implicit none
  integer, intent(in) :: ox, oy, cx, cy
  double precision, intent(in), dimension(ox,oy) :: omap
  double precision, intent(out), dimension(cx,cy) :: cmap
  integer :: nx, ny

  do nx = 1, cx
    do ny = 1, cy
      cmap(nx,ny) = omap(nx+(ox-cx)*0.5,ny+(oy-cy)*0.5)
    end do
  end do
 
end subroutine cutmap


end module utils

