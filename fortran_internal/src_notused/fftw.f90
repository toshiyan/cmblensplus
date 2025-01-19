!////////////////////////////////////////////////////!
! * dft module
!////////////////////////////////////////////////////!

module fftw
  use constants, only: dlc, iu, pi, twopi
  use utils, only: spin_weight, gaussian_alm
  implicit none

  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)

  private FFTW_ESTIMATE
  private dlc, iu, pi, twopi
  private spin_weight, gaussian_alm

contains 

!//////////////////////////////////////////////////////////////////////!
! * QU <-> EB transform in Fourier space


subroutine QU2EB(nn,D,QU,EB,trans)
! * trans = 1  :  QU -> EB
! * trans = -1 :  EB -> QU
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans
  double precision, intent(in) :: D(2)
  double complex, intent(inout), dimension(:,:) :: QU,EB
  !internal
  integer :: i, j, n
  double precision :: x, y, sin2t, cos2t

  n = 1
  do i = 1, nn(1)
    x = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if (x==0d0.and.y==0d0) then
        cos2t = 0d0
        sin2t = 0d0
      else
        cos2t = (x**2-y**2)/(x**2+y**2)
        sin2t = 2d0*x*y/(x**2+y**2)
      end if
      if (trans==1) then
        EB(1,n) = QU(n,1)*cos2t + QU(n,2)*sin2t
        EB(2,n) = -QU(n,1)*sin2t + QU(n,2)*cos2t
      else if (trans==-1) then
        QU(n,1) = EB(1,n)*cos2t - EB(2,n)*sin2t
        QU(n,2) = EB(1,n)*sin2t + EB(2,n)*cos2t
      end if
      n = n + 1
    end do
  end do 

end subroutine QU2EB


subroutine dft_all_1darray(TQU,nn,D,TEB,trans)
  !trans=-1 : T,E,B -> T,Q+iU
  !trans=+1 : T,Q+iU -> T,E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: TQU(:,:)
  complex(dlc), intent(inout) :: TEB(:,:)
  !internal
  integer :: npix
  complex(dlc), allocatable :: T(:), cP(:), P(:)

  npix = nn(1)*nn(2)

  select case(trans)
  case(-1)
    allocate(T(npix),P(npix))
    T = TEB(1,:)
    call dft(T,nn,D,-1)
    TQU(:,1) = dble(T)
    P = TEB(2,:) + iu*TEB(3,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    TQU(:,2) = dble(P)
    TQU(:,3) = aimag(P)
    deallocate(T,P)
  case(1)
    allocate(T(npix),P(npix),cP(npix))
    T  = TQU(:,1)
    call dft(T,nn,D,-1)
    TQU(:,1) = dble(T)
    P  = TQU(:,2) + iu*TQU(:,3)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call dft(cP,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call spin_weight(cP,nn,D,1,-2)
    TEB(1,:) = (P+cP)*0.5d0
    TEB(2,:) = -iu*(P-cP)*0.5d0
    deallocate(t,P,cP)
  end select

end subroutine dft_all_1darray

subroutine dft_all_2darray(TQU,nn,D,TEB,trans)
  !trans=-1 : T,E,B -> T,Q+iU
  !trans=+1 : T,Q+iU -> T,E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: TQU(:,:,:)
  complex(dlc), intent(inout) :: TEB(:,:,:)
  !internal
  complex(dlc), allocatable :: T(:,:), cP(:,:), P(:,:)

  select case(trans)
  case(-1)
    allocate(T(nn(1),nn(2)),P(nn(1),nn(2)))
    T = TEB(1,:,:)
    call dft(T,nn,D,-1)
    TQU(:,:,1) = dble(T)
    P = TEB(2,:,:) + iu*TEB(3,:,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    TQU(:,:,2) = dble(P)
    TQU(:,:,3) = aimag(P)
    deallocate(T,P)
  case(1)
    allocate(T(nn(1),nn(2)),P(nn(1),nn(2)),cP(nn(1),nn(2)))
    T  = TQU(:,:,1)
    call dft(T,nn,D,-1)
    TQU(:,:,1) = dble(T)
    P  = TQU(:,:,2) + iu*TQU(:,:,3)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call dft(cP,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call spin_weight(cP,nn,D,1,-2)
    TEB(1,:,:) = (P+cP)*0.5d0
    TEB(2,:,:) = -iu*(P-cP)*0.5d0
    deallocate(t,P,cP)
  end select

end subroutine dft_all_2darray


!///////////////////////////////////////////////////////////////////////!
!* generate Gaussian random fluctuations on 2D map

subroutine gaussian_map(nn,D,eL,Cl,map,fix)
!simulate Gaussian alm
  implicit none
  !I/O
  logical, intent(in), optional :: fix
  integer, intent(in) :: nn(1:2), eL(1:2)
  double precision, intent(in) :: D(1:2), Cl(:)
  complex(dlc), intent(out) :: map(:)
  !internal
  integer :: iL(1:2)
  complex(dlc), dimension(:), allocatable :: alm

  allocate(alm(nn(1)*nn(2)))

  if(present(fix)) then
    call Gaussian_Alm(nn,D,eL,alm,Cl,fix)
  else
    write(*,*) 'generate gaussian alm'
    call Gaussian_Alm(nn,D,eL,alm,Cl)
  end if

  write(*,*) 'alm -> map'
  call dft(alm,nn,D,-1)
  map = alm
  deallocate(alm)

end subroutine gaussian_map


subroutine gaussian_map_pol(nn,D,eL,EE,BB,QU,P,fix)
!simulate Gaussian Elm and Blm
  implicit none
  !I/O
  logical, intent(in), optional :: fix
  integer, intent(in) :: nn(1:2), eL(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(in), optional :: EE(:), BB(:)
  double precision, intent(out), optional :: QU(:,:)
  complex(dlc), intent(out), optional :: P(:)
  !internal
  double precision, dimension(:,:), allocatable :: QUfull
  complex(dlc), dimension(:,:), allocatable :: alm

  allocate(alm(2,nn(1)*nn(2)))
  alm = 0d0

  if(present(fix)) then
    if(present(EE)) call gaussian_alm(nn,D,eL,alm(1,:),EE,fix)
    if(present(BB)) call gaussian_alm(nn,D,eL,alm(2,:),BB,fix)
  else
    if(present(EE)) call gaussian_alm(nn,D,eL,alm(1,:),EE)
    if(present(BB)) call gaussian_alm(nn,D,eL,alm(2,:),BB)
  end if

  allocate(QUfull(nn(1)*nn(2),2))
  call dft_pol_1darray(QUfull,nn,D,alm,-1)
  deallocate(alm)

  if(present(QU)) QU = QUfull
  if(present(P))  P  = QUfull(:,1) + iu*QUfull(:,2)
  deallocate(QUfull)

end subroutine gaussian_map_pol


subroutine array_deriv(nn,D,W,Wx,Wy,Wxx,Wxy,Wyy)
! for apodization window
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(:) :: W
  double complex, intent(inout):: Wx(:), Wy(:), Wxx(:),Wxy(:), Wyy(:)
  integer :: i, j, n, npix
  double precision :: lx, ly
  double complex, allocatable :: Wl(:)

  npix = nn(1)*nn(2)

  !* array in Fourier space
  allocate(Wl(npix))
  Wl = W
  call dft(Wl,nn,D,1)
  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      Wx(n) = iu*lx*Wl(n)
      Wy(n) = iu*ly*Wl(n)
      Wxx(n) = -lx**2*Wl(n)
      Wxy(n) = -lx*ly*Wl(n)
      Wyy(n) = -ly**2*Wl(n)
      n = n + 1
    end do
  end do 
  deallocate(Wl)

  !* Derivatives of windows functions
  call dft(Wx,nn,D,-1)
  call dft(Wy,nn,D,-1)
  call dft(Wxx,nn,D,-1)
  call dft(Wxy,nn,D,-1)
  call dft(Wyy,nn,D,-1)

end subroutine array_deriv




end module myfftw


