module remapping
  use funcs,    only: factorial
  use myconst,  only: dlc, iu, twopi
  use anaflat,  only: rotation
  use myfftw,   only: dft, dft_pol, derivemap
  implicit none

  interface remap
    module procedure remap_1d, remap_2d
  end interface

  private factorial
  private dlc, iu, twopi
  private rotation
  private dft, dft_pol, derivemap

contains


subroutine find_correct_phi(nn,D,alpha,beta)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), alpha(:,:), beta(:,:)
  !internal
  integer :: i, n, npix
  double precision :: det
  double precision, allocatable :: bi(:,:), ai(:,:), dai(:,:,:), da(:,:,:)
  
  npix = nn(1)*nn(2)
  allocate(bi(npix,2),ai(npix,2),dai(npix,2,2),da(npix,2,2)); dai=0d0; da=0d0

  !initial
  ai = alpha
  bi = 0d0
  call derivemap(nn,D,alpha(:,1),da(:,:,1),1)
  call derivemap(nn,D,alpha(:,2),da(:,:,2),1)

  !iterate
  do i = 1, 100
    ai  = alpha
    dai = da
    call remap_1d(nn,D,bi,ai(:,1),4)
    call remap_1d(nn,D,bi,ai(:,2),4)
    call remap_1d(nn,D,bi,dai(:,1,1),4)
    call remap_1d(nn,D,bi,dai(:,1,2),4)
    call remap_1d(nn,D,bi,dai(:,2,1),4)
    call remap_1d(nn,D,bi,dai(:,2,2),4)
    do n = 1, npix
      det = (1d0-dai(n,1,1))*(1d0-dai(n,2,2)) - dai(n,1,2)*dai(n,2,1)
      bi(n,1) = bi(n,1) - ( (1d0-dai(n,2,2)) * (bi(n,1)-ai(n,1)) - dai(n,1,2) * (bi(n,2)-ai(n,2)) )
      bi(n,2) = bi(n,2) - ( (1d0-dai(n,1,1)) * (bi(n,2)-ai(n,2)) - dai(n,2,1) * (bi(n,1)-ai(n,1)) )
    end do
  end do

  deallocate(bi,ai,dai)

end subroutine find_correct_phi


subroutine rotation_alm(nn,D,alm,rot,rtype)
  implicit none
  character(*), intent(in) :: rtype
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), rot(:)
  complex(dlc), intent(inout) :: alm(:,:)
  double precision, allocatable :: QU(:,:)

  allocate(QU(nn(1)*nn(2),2))
  call dft_pol(QU,nn,D,alm,-1)
  call rotation(QU,rot,rtype)
  call dft_pol(QU,nn,D,alm,1)
  deallocate(QU)

end subroutine rotation_alm


subroutine kap2def(nn,D,kap,def,lcut)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(in), optional :: lcut
  complex(dlc), intent(in) :: kap(:,:)
  complex(dlc), intent(out) :: def(:,:,:)
  !internal
  integer :: a, b
  double precision :: lx, ly, lc

  lc = 8000d0
  if (present(lcut)) lc = lcut

  do a = 1, nn(1)
    do b = 1, nn(2)
      lx = twopi*dble(a-1-nn(1)*0.5d0)/D(1)
      ly = twopi*dble(b-1-nn(2)*0.5d0)/D(2)
      if (lx==0d0.and.ly==0d0) cycle
      if (dsqrt(lx**2+ly**2)>lc) cycle !multipole cut
      def(a,b,1) = 2d0*iu*lx*kap(a,b)/(lx**2+ly**2) ! (-2/L^2) x (-iu L) kappa
      def(a,b,2) = 2d0*iu*ly*kap(a,b)/(lx**2+ly**2) 
    end do
  end do
  call dft(def(:,:,1),nn,D,-1)
  call dft(def(:,:,2),nn,D,-1)

end subroutine kap2def


subroutine kap2def_1d(nn,D,kap,def,lcut)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(in), optional :: lcut
  complex(dlc), intent(in) :: kap(:)
  complex(dlc), intent(out) :: def(:,:)
  !internal
  integer :: a, b, n
  double precision :: lx, ly, lc

  lc = 8000d0
  if (present(lcut)) lc = lcut

  n = 0
  do a = 1, nn(1)
    do b = 1, nn(2)
      lx = twopi*dble(a-1-nn(1)*0.5d0)/D(1)
      ly = twopi*dble(b-1-nn(2)*0.5d0)/D(2)
      n = n + 1
      if (lx==0d0.and.ly==0d0) cycle
      if (dsqrt(lx**2+ly**2)>lc) cycle !multipole cut
      def(n,1) = 2d0*iu*lx*kap(n)/(lx**2+ly**2) ! (-2/L^2) x (-iu L) kappa
      def(n,2) = 2d0*iu*ly*kap(n)/(lx**2+ly**2) 
    end do
  end do
  call dft(def(:,1),nn,D,-1)
  call dft(def(:,2),nn,D,-1)

end subroutine kap2def_1d


subroutine remap_lin(nn,D,lx,ly,Tlm,plm)
! remapping with only dphi x dT term
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), lx(:), ly(:)
  complex(dlc), intent(in) :: plm(:) !lensing potential
  complex(dlc), intent(inout) :: Tlm(:) !temperature Fourier mode
  !internal
  integer :: npix
  complex(dlc), allocatable :: def(:,:), dT(:,:)

  npix = nn(1)*nn(2)

  allocate(def(npix,2),dT(npix,2))

  !deflection angle
  def(:,1) = -iu*lx*plm
  def(:,2) = -iu*ly*plm
  !derivative of temperature
  dT(:,1) = -iu*lx*Tlm
  dT(:,2) = -iu*ly*Tlm

  call dft(def(:,1),nn,D,-1)
  call dft(def(:,2),nn,D,-1)
  call dft(dT(:,1),nn,D,-1)
  call dft(dT(:,2),nn,D,-1)

  !* remapping by deflection vector
  dT = def*dT
  call dft(dT(:,1),nn,D,1)
  call dft(dT(:,2),nn,D,1)
  Tlm = Tlm + sum(dT,dim=2)

  deallocate(def,dT)

end subroutine remap_lin


subroutine remap_1d(nn,D,vec,map,nth)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nth
  double precision, intent(in) :: D(1:2), vec(:,:)
  double precision, intent(inout) :: map(:)
  !internal 
  double precision, allocatable :: m2d(:,:), v2d(:,:,:)

  allocate(m2d(nn(1),nn(2)),v2d(nn(1),nn(2),2))
  m2d        = reshape(map,nn,order=[2,1])
  v2d(:,:,1) = reshape(vec(:,1),nn,order=[2,1])
  v2d(:,:,2) = reshape(vec(:,2),nn,order=[2,1])
  call remap_2d(D,v2d,m2d,nth)
  map = reshape(transpose(m2d),[nn(1)*nn(2)])
  deallocate(m2d,v2d)

end subroutine remap_1d


subroutine remap_2d(D,vec,map,nth)
  implicit none
  !I/O
  double precision, intent(in) :: D(1:2), vec(:,:,:)
  double precision, intent(inout) :: map(:,:)
  !internal
  integer :: i,j,a,b,k,n,nth,i0,j0,nn(2)
  double precision :: dx,dy,x0,y0,da(2),kk,nk,pm
  double precision, allocatable :: dmap(:,:,:,:), xi(:), yi(:)

  write(*,*) 'remap'

  !* get information
  nn(1) = size(map,dim=1)
  nn(2) = size(map,dim=2)

  !* define coordinate
  allocate(xi(nn(1)),yi(nn(2)))
  dx = D(1)/dble(nn(1))
  dy = D(2)/dble(nn(2))
  xi = [ ( (dble(i-1)), i=1,nn(1) ) ] * dx
  yi = [ ( (dble(i-1)), i=1,nn(2) ) ] * dy

  !* interpolation on regular grid
  allocate(dmap(nth+1,nth+1,nn(1),nn(2))); dmap=0d0

  !* prepare derivative maps
  write(*,*) 'higher-order derivatives'
  call derivemap(D,map,dmap)

  !* compute lensed map
  do a=1, nn(1)
    do b=1, nn(2)

      !shifted position
      x0 = xi(a)+vec(a,b,1)
      y0 = yi(b)+vec(a,b,2)

      !find nearest grid index
      pm = 1d0
      if (vec(a,b,1)<0d0) pm = -1d0
      i0 = a + int(vec(a,b,1)/dx+0.5*pm)
      pm = 1d0
      if (vec(a,b,2)<0d0) pm = -1d0
      j0 = b + int(vec(a,b,2)/dy+0.5*pm)
      if (i0<1) i0=1
      if (j0<1) j0=1
      if (i0>nn(1)) i0=nn(1)
      if (j0>nn(2)) j0=nn(2)

      !delta alpha
      da(1) = x0 - xi(i0)
      da(2) = y0 - yi(j0)

      !lensed map
      map(a,b) = dmap(1,1,i0,j0) !0th order
      do n = 1, nth !higher order
        do k = 0, n
          kk = factorial(k)
          nk = factorial(n-k)
          map(a,b) = map(a,b) + (da(1)**(n-k)*da(2)**k)/(kk*nk) * dmap(k+1,n+1,i0,j0)
        end do
      end do

    end do
  end do

  deallocate(dmap,xi,yi)

end subroutine remap_2d


!////////// (d phi)^n x d^n scal  //////////!
subroutine dphidmap_1d(nn,D,kap,alm,nth)
  implicit none
  !I/O
  integer, intent(in) :: nth, nn(2)
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in) :: kap(:)
  complex(dlc), intent(inout) :: alm(:)
  !internal
  integer :: k, npix
  double precision :: kk, nk
  double precision, allocatable :: dmap(:,:)
  complex(dlc), allocatable :: blm(:), da(:,:)

  !* get information
  npix = nn(1)*nn(2)

  !* prepare derivative maps
  allocate(dmap(nth+1,npix),blm(npix),da(npix,2)); dmap=0d0; blm=0d0; da=0d0
  call derivemap(nn,D,alm,dmap,nth)
  call kap2def_1d(nn,D,kap,da)
  do k = 0, nth
    kk  = factorial(k)
    nk  = factorial(nth-k)
    blm = blm + (da(:,1)**(nth-k)*da(:,2)**k)/(kk*nk) * dmap(k+1,:)
  end do
  call dft(blm,nn,D,1)
  alm = blm
  deallocate(dmap,blm,da)

end subroutine dphidmap_1d


subroutine dphidpol_1d(nn,D,kap,alm,nth)
  implicit none
  !I/O
  integer, intent(in) :: nth, nn(2)
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in) :: kap(:)
  complex(dlc), intent(inout) :: alm(:,:)
  !internal
  integer :: k, npix
  double precision :: kk, nk
  double precision, allocatable :: QU(:,:), dmap(:,:,:)
  complex(dlc), allocatable :: da(:,:)

  !* get information
  npix = nn(1)*nn(2)

  !* prepare derivative maps
  allocate(dmap(2,nth+1,npix),da(npix,2),QU(npix,2)); dmap=0d0; da=0d0; QU=0d0
  call dft_pol(QU,nn,D,alm,-1)
  call derivemap(nn,D,QU(:,1),dmap(1,:,:),nth)
  call derivemap(nn,D,QU(:,2),dmap(2,:,:),nth)
  call kap2def_1d(nn,D,kap,da)
  QU = 0d0
  do k = 0, nth
    kk = factorial(k)
    nk = factorial(nth-k)
    QU(:,1) = QU(:,1) + (da(:,1)**(nth-k)*da(:,2)**k)/(kk*nk) * dmap(1,k+1,:)
    QU(:,2) = QU(:,2) + (da(:,1)**(nth-k)*da(:,2)**k)/(kk*nk) * dmap(2,k+1,:)
  end do
  call dft_pol(QU,nn,D,alm,1)
  deallocate(dmap,da,QU)

end subroutine dphidpol_1d


subroutine dphi2d2T(nn,D,pT,k1,k2,d1,d2,alm,ddT)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in), optional :: k1(:), k2(:), d1(:,:), d2(:,:), alm(:)
  double precision, intent(in), optional :: ddT(:,:)
  complex(dlc), intent(out) :: pT(:)
  !internal
  integer :: k, npix
  double precision, allocatable :: dmap(:,:)
  complex(dlc), allocatable :: blm(:), da1(:,:), da2(:,:)

  !* get information
  npix = nn(1)*nn(2)

  allocate(dmap(3,npix),blm(npix),da1(npix,2),da2(npix,2)); dmap=0d0; blm=0d0; da1=0d0; da2=0d0

  !* deflection angle maps
  if (present(d1)) da1 = d1
  if (present(k1)) call kap2def_1d(nn,D,k1,da1)
  da2 = da1
  if (present(d2)) da2 = d2
  if (present(k2)) call kap2def_1d(nn,D,k2,da2)

  !* derivative T maps
  if (present(alm)) call derivemap(nn,D,alm,dmap,2)
  if (present(ddT)) dmap = ddT

  !* compute dphi^2 x d^2T in map space
  blm = ( da1(:,1)*da2(:,1)*dmap(1,:) + (da1(:,1)*da2(:,2)+da1(:,2)*da2(:,1))*dmap(2,:) + da1(:,2)*da2(:,2)*dmap(3,:) )/2d0

  !* map to alm
  call dft(blm,nn,D,1)
  pT = blm

  deallocate(dmap,blm,da1,da2)

end subroutine dphi2d2T


end module remapping

