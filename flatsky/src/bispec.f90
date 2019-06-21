!////////////////////////////////////////////////////!
! bispectrum measurement from flatsky simulations
!////////////////////////////////////////////////////!

module bispflat
  use myconst, only: dlc
  use anaflat, only: elarray
  use myfftw,  only: dft

  private dlc
  private elarray
  private dft

contains

subroutine binfilter(nn,D,bp,bf)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bp(:), D(2)
  double precision, intent(out) :: bf(:,:)
  integer :: n, b, npix, bmax
  double precision, allocatable :: els(:)

  npix = nn(1)*nn(2)
  bmax = size(bp)-1

  bf = 0d0

  allocate(els(npix))
  els = elarray(nn,D)

  do b = 1, bmax
    do n = 1, npix
      if (bp(b)<=els(n).and.els(n)<bp(b+1)) bf(b,n) = 1d0
    end do
  end do

  deallocate(els)

end subroutine binfilter


subroutine norm_fold(nn,D,bf,dlb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:,:), D(2)
  double precision, intent(out) :: dlb(:)
  integer :: b, npix
  complex(dlc), allocatable :: wlm(:,:)

  npix = nn(1)*nn(2)

  do b = 1, size(bf,dim=2)
    allocate(wlm(2,npix));  wlm=0d0
    write(*,*) b
    wlm(1,:) = bf(1,b,:)
    wlm(2,:) = bf(2,b,:)
    call dft(wlm(1,:),nn,D,-1)
    call dft(wlm(2,:),nn,D,-1)
    dlb(b) = D(1)*D(2)*sum(wlm(1,:)*wlm(2,:)**2)
    deallocate(wlm)
  end do

end subroutine norm_fold


subroutine bisp_fold(nn,D,bf,dlb,klm,Bb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:,:), D(2), dlb(:)
  double precision, intent(out) :: Bb(:)
  complex(dlc), intent(in) :: klm(:)
  integer :: b, npix
  complex(dlc), allocatable :: wklm(:,:)

  npix = nn(1)*nn(2)
  
  do b = 1, size(bf,dim=2)
    allocate(wklm(2,npix));  wklm=0d0
    write(*,*) 'binned kappa', b
    wklm(1,:) = bf(1,b,:)*klm
    wklm(2,:) = bf(2,b,:)*klm
    call dft(wklm(1,:),nn,D,-1)
    call dft(wklm(2,:),nn,D,-1)
    Bb(b) = sum(wklm(1,:)*wklm(2,:)**2)/dlb(b)
    deallocate(wklm)
  end do

end subroutine bisp_fold


subroutine norm_equi(nn,D,bf,dlb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:), D(2)
  double precision, intent(out) :: dlb(:)
  integer :: b, npix
  complex(dlc), allocatable :: wlm(:)

  npix = nn(1)*nn(2)

  do b = 1, size(bf,dim=1)
    allocate(wlm(npix));  wlm=0d0
    write(*,*) b
    wlm = bf(b,:)
    call dft(wlm,nn,D,-1)
    !dlb(1,b) = D(1)*D(2)*sum(wlm**2) !power spectrum
    dlb(b) = D(1)*D(2)*sum(wlm**3) !bispectrum 
    deallocate(wlm)
  end do

end subroutine norm_equi


subroutine bisp_equi(nn,D,bf,dlb,klm,Bb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:), D(2), dlb(:)
  double precision, intent(out) :: Bb(:)
  complex(dlc), intent(in) :: klm(:)
  integer :: b, npix
  complex(dlc), allocatable :: wklm(:)

  npix = nn(1)*nn(2)
  
  do b = 1, size(bf,dim=1)
    allocate(wklm(npix));  wklm=0d0
    write(*,*) 'binned kappa', b
    wklm = bf(b,:)*klm
    call dft(wklm,nn,D,-1)
    !Bb(1,b) = sum(wklm**2)/dlb(1,b)
    Bb(b) = sum(wklm**3)/dlb(b)
    deallocate(wklm)
  end do

end subroutine bisp_equi


end module bispflat

