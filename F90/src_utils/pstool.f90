!////////////////////////////////////////////////////!
! * Angular power spectrum calculation
!////////////////////////////////////////////////////!

module pstool
  use constants, only: pi, TT, TE, EE, BB, dd, Td, Ed, oo
  use general,   only: splint, spline, FileColumns, FileLines, savetxt, linspace, interp_lin, check_positive
  implicit none

  interface power_binning
    module procedure :: power_binning_full, power_binning_flat_1d, power_binning_flat_2d
  end interface

  private pi, TT, TE, EE, BB, dd, Td, Ed, oo
  private splint, spline, FileColumns, FileLines, savetxt, linspace, interp_lin, check_positive

contains


!//////////////////////////////////////////////////////////////////////////////!
! Binned cls
!//////////////////////////////////////////////////////////////////////////////!

subroutine binned_ells(eL,bp,bc,spc)
! return binned multipole edges and centers
  implicit none
  !I/O
  character(*), intent(in), optional :: spc
  integer, intent(in) :: eL(:)
  double precision, intent(out), optional :: bp(:), bc(:)
  !internal
  character(8) :: sp=''
  integer :: i, n
  double precision :: d, edge(2)
  double precision, allocatable :: p(:)

  if (eL(1)<=0.and.sp/='') stop 'ell minimum should be > 0 for spacing'

  if (present(bc))   n = size(bc)
  if (present(bp))   n = size(bp) - 1
  if (present(spc))  sp = spc

  edge = dble(eL)

  allocate(p(0:n))
  p(0) = eL(1)

  select case(sp)
  case('log')
    d     = dlog(edge(2)/edge(1))/dble(n)
    p(1:) = [ ( ( edge(1)*dexp(d*i) ), i=1,n ) ]
  case('log10')
    d     = dlog10(edge(2)/edge(1))/dble(n)
    p(1:) = [ ( ( edge(1)*10**(d*i) ), i=1,n ) ]
  case('p2')
    d     = (dsqrt(edge(2))-dsqrt(edge(1)))/dble(n)
    p(1:) = [ ( ( (dsqrt(edge(1))+d*i)**2) ,i=1,n ) ]
  case('p3')
    d     = (edge(2)**(1d0/3d0)-edge(1)**(1d0/3d0))/dble(n)
    p(1:) = [ ( ( (edge(1)**3+d*i)**3 ), i=1,n ) ]
  case default
    d     = (edge(2)-edge(1))/dble(n)
    p(1:) = [ ( ( edge(1)+d*i ), i=1,n ) ]
  end select

  if (present(bc)) bc = [(((p(i)+p(i-1))*0.5d0),i=1,n)]
  if (present(bp)) bp = p
 
  deallocate(p)


end subroutine binned_ells


!//////////////////////////////////////////////////////////////////////////////!
! interpretation of binned cls to unbinned cls
!//////////////////////////////////////////////////////////////////////////////!

subroutine cb2cl(bc,Cb,Cl,bp,method)
! interpolate non-binned Cl from binned Cl
  implicit none
  !I/O
  double precision, intent(in) :: bc(:), Cb(:)
  double precision, intent(out) :: Cl(:)
  !(optional)
  character(*), intent(in), optional :: method
  double precision, intent(in), optional :: bp(:)
  !internal
  logical :: p
  character(8) :: m
  integer :: ln, bn, i, j
  double precision :: y2(size(Cb)), Cb0(size(Cb)+1)

  m  = ''
  bn = size(Cb)
  Cl = 0d0
  ln = size(Cl)
  if (present(method)) m = method

  call check_positive(Cb,p) !check positivity

  !interpolation
  select case (m)
  case ('step')
    if (.not.present(bp)) stop '(cb2cl): bin edge is required'
    do i = 1, bn
      if (ln<int(bp(i))) goto 1
      do j = int(bp(i)), min(ln,int(bp(i+1)))
        Cl(j) = Cb(i)
      end do
    end do
  case ('linear')
    !edge 1
    !Cl(1:int(bc(1))) = Cb(1)
    do j = 1, int(bc(1))
      Cl(j) = max(0d0,interp_lin(dble(j),bc(1),bc(2),Cb(1),Cb(2)))
    end do
    do i = 1, bn-1
      if (ln<int(bc(i))) goto 1
      do j = int(bc(i)), min(ln,int(bc(i+1)))
        Cl(j) = interp_lin(dble(j),bc(i),bc(i+1),Cb(i),Cb(i+1))
      end do
    end do
    !edge 2
    !Cl(int(bc(bn))+1:) = Cb(bn)
    do j = int(bc(bn))+1, ln
      Cl(j) = max(0d0,interp_lin(dble(j),bc(bn-1),bc(bn),Cb(bn-1),Cb(bn)))
    end do
  case default
    call spline(bc,Cb,bn,0d0,0d0,y2)
    if (ln<int(bc(1))) goto 1
    do j = int(bc(1)), min(ln,int(bc(bn)))
      Cl(j) = splint(dble(j),bc,Cb,y2)
    end do
1 end select

  !check cl
  do j = 1, size(Cl)
    if (Cl(j)>=0d0.or..not.p) cycle
    write(*,*) 'warning(cb2cl): interpolated power is negative', j, Cl(j)
    Cl(j) = 0d0
  end do

end subroutine cb2cl


subroutine cl_interp_spline(bls,bCl,tL,Cl,islog)
  implicit none
  !I/O
  logical, intent(in) :: islog 
  integer, intent(in) :: tL(:)
  double precision, intent(in) :: bls(:), bCl(:)
  double precision, intent(out) :: Cl(:)
  !internal
  integer :: l, n, i, nmax
  double precision, allocatable :: tbCl(:), tbls(:), tbCldd(:)

  !choose non-zero bCls
  n = 0
  do i = 1, size(bCl)
    if(bCl(i)>0) n = n + 1
  end do
  nmax = n

  allocate(tbCl(nmax),tbls(nmax),tbCldd(nmax))
  n = 0
  do i = 1, size(bCl)
    if(bCl(i)<=0) cycle
    n = n + 1
    tbCl(n) = bCl(i)
    tbls(n) = bls(i)
  end do

  Cl = 0d0
  call spline(tbls,tbCl,nmax,0d0,0d0,tbCldd)
  do l = tL(1), tL(2)
    if(islog) then 
      Cl(l) = splint(dlog(dble(l)),tbls,tbCl,tbCldd)
    else 
      Cl(l) = splint(dble(l),tbls,tbCl,tbCldd)
    end if
  end do

  deallocate(tbCl,tbls)

end subroutine cl_interp_spline


!//////////////////////////////////////////////////////////////////////////////!
! Binning angular power spectrum
!//////////////////////////////////////////////////////////////////////////////!


subroutine power_binning_flat_1d(bp,els,AL,BL,bA,bV)
!* compute averaging factor: 
!   bA = bV^2 * \sum_{l=bp(i)}^{bp(i+1)} BL(l)*AL(l)
! AL: amplitude
! BL: band-pass filter
!
  implicit none
!
! [input]
!   bp  --- edges of multipole binning
!   els --- multipoles as a function of npix
!   AL  --- amplitude
!   BL  --- band-pass filter
  double precision, intent(in) :: bp(:), els(:), AL(:), BL(:)
!
! [output]
!   bA  --- binned amplitude
!   bV  --- binned variance
  double precision, intent(out) :: bA(:), bV(:)
!
! [internal]
  integer :: b, i, npix, bmax
  double precision, allocatable :: bW(:)

  npix = size(els)
  bmax = size(bp) - 1

  allocate(bW(bmax))
  !* bW = [ sum_i BL(i) ]^-1 = bV^2 : optimal variance at b-th bin
  bW = 0d0;  bA = 0d0;  bV = 0d0
  do b = 1, bmax
    do i = 1, npix
      if (els(i)<bp(b).or.els(i)>=bp(b+1)) cycle
      bW(b) = bW(b) + BL(i) 
      bA(b) = bA(b) + BL(i)*Al(i)
    end do
    !if (bW(b) <= 0d0) write(*,*) 'no grid data points are assigned to this ell-bin ...'
    if (bW(b) <= 0d0) cycle
    bW(b) = 1d0/bW(b) !square of variance
    bA(b) = bA(b)*bW(b)
    bV(b) = dsqrt(bW(b))
  end do

  deallocate(bW)

end subroutine power_binning_flat_1d


subroutine power_binning_flat_2d(bp,els,AL,BL,bA,bV)
!* compute averaging factor: 
!   bA = bV^2 * \sum_{l=bp(i)}^{bp(i+1)} BL(l)*AL(l)
! AL: amplitude
! BL: band-pass filter
!
  implicit none
!
! [input]
!   bp  --- edges of multipole binning
!   els --- multipoles as a function of npix
!   AL  --- amplitude
!   BL  --- band-pass filter
  double precision, intent(in) :: bp(:), els(:,:), AL(:,:), BL(:,:)
!
! [output]
!   bA  --- binned amplitude
!   bV  --- binned variance
  double precision, intent(out) :: bA(:), bV(:)
!
! [internal]
  integer :: b, i, j, nx, ny, bmax
  double precision, allocatable :: bW(:)

  nx   = size(els,dim=1)
  ny   = size(els,dim=2)
  bmax = size(bp) - 1

  allocate(bW(bmax))
  !* bW = [ sum_i BL(i) ]^-1 = bV^2 : optimal variance at b-th bin
  bW = 0d0;  bA = 0d0;  bV = 0d0
  do b = 1, bmax
    do i = 1, nx
      do j = 1, ny
        if (els(i,j)<bp(b).or.els(i,j)>=bp(b+1)) cycle
        bW(b) = bW(b) + BL(i,j) 
        bA(b) = bA(b) + BL(i,j)*Al(i,j)
      end do
    end do
    !if (bW(b) <= 0d0) write(*,*) 'no grid data points are assigned to this ell-bin ...'
    if (bW(b) <= 0d0) cycle
    bW(b) = 1d0/bW(b) !square of variance
    bA(b) = bA(b)*bW(b)
    bV(b) = dsqrt(bW(b))
  end do

  deallocate(bW)

end subroutine power_binning_flat_2d


!//////////////////////////////////////////////////////////////////////////////!
! From alm to angular power spectrum (fullsky)
!//////////////////////////////////////////////////////////////////////////////!

subroutine power_binning_full(bp,eL,Cl,bC,Vl)
  implicit none 
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: bp(:), Cl(:)
  double precision, intent(in), optional :: Vl(:)
  double precision, intent(out) :: bC(:)
  !internal 
  integer :: b, i, bmax
  double precision, allocatable :: bW(:), Wl(:)

  bmax = size(bp) - 1

  allocate(bW(bmax),Wl(size(Cl))); bW=0d0; Wl=1d0
  !* compute optimal variance of b-th bin, bW(b)
  if (present(Vl))  Wl = Vl
  bC = 0d0
  do b = 1, bmax
    do i = el(1), el(2)
      if (i<bp(b).or.i>=bp(b+1)) cycle
      bW(b) = bW(b) + Wl(i)
      bC(b) = bC(b) + Wl(i)*Cl(i)
    end do
    if (bW(b)<=0d0) cycle
    bW(b) = 1d0/bW(b)
    bC(b) = bC(b)*bW(b)
  end do
  deallocate(bW,Wl)

end subroutine power_binning_full


subroutine readcl_camb(cl,f,eL,bb,lsq)
!* Read cls from CAMB output files. This routine assumes that the file is obtained by CAMB and contains 
!*
!*   bb = False: l, TT, EE, TE, dd, Td, (Ed)
!*   bb = True:  l, TT, EE, BB, TE, (curl)
!*
!* lsq = True to remove multipole factor
!*
  implicit none 
![input]
  character(*), intent(in) :: f
  integer, intent(in) :: eL(1:2)
  double precision, intent(out) :: Cl(:,:)
  logical, intent(in) :: lsq, bb
![internal]
  double precision :: rC(1:10)
  integer :: i, l, n, m, k, nn

  k = size(cl,dim=1)

  open(unit=20,file=trim(f),status='old')

  n = FileColumns(20)
  m = FileLines(20)

  write(*,*) n, m

  do i = 1, m

    read(20,*) l, rC(1:n-1)

    if (el(1)>l.or.l>el(2)) cycle

    if (lsq) then
      rC = rC*2d0*pi/dble(l**2+l)
      if (.not.bb) then
        if(n>=5) rC(4) = rC(4)*dble(l+1d0)/(2d0*pi*dble(l)**3)
        if(n>=6) rC(5) = rC(5)*dble(l+1d0)/(2d0*pi*dble(l)**2)
      end if
    end if

    ! to output
    do nn = 2, n
      if(k>=nn-1) cl(nn-1,l) = rC(nn-1)
    end do

  end do

  close(20)

end subroutine readcl_camb


subroutine map_vars(sigma,eL,cl)
! * variance of a map and its derivative
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: cl(:)
  double precision, intent(out) :: sigma(0:1)
  !internal
  integer :: l
  double precision :: al
  
  sigma = 0d0
  do l = el(1), el(2)
    al = dble(l)
    sigma(0) = sigma(0) + (2d0*al+1d0)*cl(l)
    sigma(1) = sigma(1) + (2d0*al+1d0)*(al**2+al)*cl(l)
  end do 
  sigma = sigma/(4d0*pi)

end subroutine map_vars


end module pstool

