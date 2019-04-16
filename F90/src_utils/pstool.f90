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
! * return binned multipole edges and centers
  implicit none
  !I/O
  character(*), intent(in), optional :: spc
  integer, intent(in) :: eL(:)
  double precision, intent(out), optional :: bp(:), bc(:)
  !internal
  character(8) :: spcs=''
  integer :: i, n
  double precision :: d
  double precision, allocatable :: p(:)

  if (present(bc))  n = size(bc)
  if (present(bp))  n = size(bp) - 1
  if (present(spc)) spcs = spc

  allocate(p(n+1))
  p(1) = eL(1)

  select case(spcs)
  case('log')
    d     = dlog(dble(eL(2))/dble(eL(1)))/dble(n)
    p(2:) = [((eL(1)*exp(d*(i-1))),i=2,n+1)]
  case default
    d     = dble(eL(2)-eL(1))/dble(n)
    p(2:) = [((d*(i-1)+p(1)),i=2,n+1)]
  end select
  if (present(bc)) bc = [(((p(i+1)+p(i))*0.5d0),i=1,n)]
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


subroutine cl2cb(bc,cl,cb)
  implicit none
  !I/O
  double precision, intent(in) :: cl(:), bc(:)
  double precision, intent(out) :: cb(:)
  !internal
  integer :: i

  cb = 0d0
  do i = 1, size(bc)
    cb(i) = cl(int(bc(i)))
  end do

end subroutine cl2cb


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


subroutine readcl_camb(Cl,f,eL,HasBB,rawCl)
!* Read cls from CAMB output files
!* This routine assumes that the file is obtained by CAMB and contains 
!*   l, TT, EE, TE, dd, Td, (Ed) 
!* if hasbb = false (e.g., scalcls.dat), or 
!*   l, TT, EE, BB, TE, (curl)
!* if hasbb = true (e.g., lensedcls.dat). 
!* rawCl is for multipole factor
  implicit none 
!
![input]
  character(*), intent(in) :: f
  integer, intent(in) :: eL(1:2)
  double precision, intent(out) :: Cl(:,:)
!
!(optional)
  logical, intent(in), optional :: hasbb, rawCl
!
![internal]
  double precision :: rC(1:10)
  integer :: i, l, n, m

  open(unit=20,file=trim(f),status='old')

  n = FileColumns(20)
  m = FileLines(20)

  do i = 1, m

    read(20,*) l, rC(1:n-1)

    if (el(1)>l.or.l>el(2)) cycle

    if (.not.(present(rawCl).and.rawCl)) then
      rC = rC*2d0*pi/dble(l**2+l)
      if (.not.(present(hasbb).and.hasbb)) then
        rC(4) = rC(4)*dble(l+1d0)/(2d0*pi*dble(l)**3)
        rC(5) = rC(5)*dble(l+1d0)/(2d0*pi*dble(l)**2)
      end if
    end if

    cl(1,l) = rC(1) !TT
    cl(2,l) = rC(2) !EE

    if(present(hasbb).and.hasbb) then
      cl(3,l) = rC(3) !BB
      cl(4,l) = rC(4) !TE
    else
      cl(3,l) = rC(3) !TE
      cl(4,l) = rC(4) !dd
      cl(5,l) = rC(5) !Td
      if(n>=7) cl(6,l) = rC(6) !Ed
    end if

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

