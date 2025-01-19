!///////////////////////////////////////////////////////////!
! * Delensing kernels
!///////////////////////////////////////////////////////////!

module delens
  use constants, only: iu, i1, pi, dlc, twopi
  use utils, only: elarrays
  use ffttools,  only: dft
  implicit none

  interface lensingb
    module procedure lensingb_from_eb, lensingb_from_e
  end interface lensingb

  private iu, i1, pi, dlc, twopi
  private elarrays
  private dft

contains 


subroutine LensingT(nn,D,T,glm,rL,Tl,gtype)
!* Return dT*dphi (-int dl l1*l2 Tl1 pl2)
  implicit none
!
! [input]
!   nn  --- x and y grids
!   rL  --- multipole ranges for T and phi
!   D   --- map size
!   T   --- temperature
!   glm --- phi map
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in) :: T(:), glm(:)
!
! (optional)
  character(*), intent(in), optional :: gtype
!
! [output]
!   Tl --- estimated lensing temperature
  complex(dlc), intent(out) :: Tl(:)
!
! [internal]
  integer :: n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: phi(:,:), alm(:), bT(:,:)

  npix = nn(1)*nn(2)
  allocate(phi(2,npix),bT(2,npix),l(2,npix),els(npix),li(npix));  phi=0d0;  bT=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),eli=li,els=els)

  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    phi(1:2,n) = l(1:2,n)*glm(n)
    bT(1:2,n)  = l(1:2,n)*T(n)
  end do 
  if(present(gtype).and.gtype=='k')  then
    phi(1,:) = phi(1,:)*(-2d0)*li**2
    phi(2,:) = phi(2,:)*(-2d0)*li**2
  end if
  deallocate(li)

  !* convolution
  call dft(phi(1,:),nn,D,-1)
  call dft(phi(2,:),nn,D,-1)
  call dft(bT(1,:),nn,D,-1)
  call dft(bT(2,:),nn,D,-1)
  allocate(alm(npix))
  alm = phi(1,:)*bT(1,:) + phi(2,:)*bT(2,:)
  deallocate(phi,bT)
  call dft(alm,nn,D,1)
  Tl = -alm
  deallocate(alm)

end subroutine LensingT


subroutine lensingb_from_eb(nn,D,alm,glm,rL,lblm,gtype)
!* Return lensing B mode using unlensed E and B modes
  implicit none
!
! [input]
!   nn  --- x and y grids
!   rL  --- multipole ranges for E and phi
!   D   --- map size
!   alm --- E and B modes
!   glm --- phi map
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in) :: alm(:,:), glm(:)
!
! (optional)
  character(*), intent(in), optional :: gtype
!
! [output]
!   lblm --- estimated lensing B-mode
  complex(dlc), intent(out) :: lblm(:)
!
! [internal]
  integer :: n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: phi(:,:), zlm(:,:), bE(:,:), bB(:,:), ei2p(:)

  npix = nn(1)*nn(2)

  allocate(phi(2,npix),bE(4,npix),bB(4,npix),ei2p(npix),l(2,npix),els(npix),li(npix))
  phi=0d0;  bE=0d0;  bB=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,eli=li,ei2p=ei2p)

  do n = 1, npix
    if (rL(1)>els(n).or.els(n)>rL(2)) cycle
    phi(1:2,n) = l(1:2,n)*glm(n)
    bE(1:2,n)  = l(1:2,n)*alm(1,n)*ei2p(n)
    bE(3:4,n)  = l(1:2,n)*alm(1,n)*conjg(ei2p(n))
    bB(1:2,n)  = l(1:2,n)*alm(2,n)*ei2p(n)
    bB(3:4,n)  = l(1:2,n)*alm(2,n)*conjg(ei2p(n))
  end do 
  if(present(gtype).and.gtype=='k')  then
    phi(1,:) = phi(1,:)*(-2d0)*li**2
    phi(2,:) = phi(2,:)*(-2d0)*li**2
  end if

  !* convolution
  call dft(phi(1,:),nn,D,-1)
  call dft(phi(2,:),nn,D,-1)
  do n = 1, 4
    call dft(bE(n,:),nn,D,-1)
    call dft(bB(n,:),nn,D,-1)
  end do
  allocate(zlm(4,npix))
  zlm(1,:) = phi(1,:)*bE(1,:) + phi(2,:)*bE(2,:)
  zlm(2,:) = phi(1,:)*bE(3,:) + phi(2,:)*bE(4,:)
  zlm(3,:) = phi(1,:)*bB(1,:) + phi(2,:)*bB(2,:)
  zlm(4,:) = phi(1,:)*bB(3,:) + phi(2,:)*bB(4,:)
  deallocate(phi,bE,bB)
  do n= 1, 4
    call dft(zlm(n,:),nn,D,1)
  end do
  lblm = -(zlm(1,:)*conjg(ei2p)-zlm(2,:)*ei2p)/(2d0*iu) - (zlm(3,:)*conjg(ei2p)+zlm(4,:)*ei2p)/2d0
  deallocate(zlm,ei2p,l,els,li)

end subroutine lensingb_from_eb


subroutine lensingb_from_e(nn,D,elm,glm,rL,lblm,gtype)
!* Return lensing B mode from unlensed E mode alone
  implicit none
!
! [input]
!   nn  --- x and y grids
!   rL  --- multipole ranges for E and phi
!   D   --- map size
!   elm --- E mode
!   glm --- phi map
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in) :: elm(:), glm(:)
!
! (optional)
  character(*), intent(in), optional :: gtype
!
! [output]
!   lblm --- estimated lensing B mode
  complex(dlc), intent(out) :: lblm(:)
!
! [internal]
  integer :: n, npix
  double precision, allocatable :: l(:,:), els(:), li(:)
  complex(dlc), allocatable :: phi(:,:), alm(:,:), bE(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(phi(2,npix),be(4,npix),ei2p(npix),l(2,npix),els(npix),li(npix));  phi=0d0;  bE=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,eli=li,ei2p=ei2p)

  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    phi(1:2,n) = l(1:2,n)*glm(n)
    bE(1:2,n)  = l(1:2,n)*elm(n)*ei2p(n)
    bE(3:4,n)  = l(1:2,n)*elm(n)*conjg(ei2p(n))
  end do 
  if(present(gtype).and.gtype=='k')  then
    phi(1,:) = phi(1,:)*(-2d0)*li**2
    phi(2,:) = phi(2,:)*(-2d0)*li**2
  end if

  !* convolution
  call dft(phi(1,:),nn,D,-1)
  call dft(phi(2,:),nn,D,-1)
  do n = 1, 4
    call dft(bE(n,:),nn,D,-1)
  end do
  allocate(alm(2,npix))
  alm(1,:) = phi(1,:)*bE(1,:) + phi(2,:)*bE(2,:)
  alm(2,:) = phi(1,:)*bE(3,:) + phi(2,:)*bE(4,:)
  deallocate(phi,bE)
  call dft(alm(1,:),nn,D,1)
  call dft(alm(2,:),nn,D,1)
  lblm = -(alm(1,:)*conjg(ei2p)-alm(2,:)*ei2p)/(2d0*iu)
  deallocate(alm,ei2p,l,els,li)

end subroutine lensingb_from_e


end module delens


