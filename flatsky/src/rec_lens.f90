!///////////////////////////////////////////////////////////!
! * Lensing Reconstruction Kernel
!///////////////////////////////////////////////////////////!

module rec_lens
  use constants, only: iu, pi, twopi
  use grid2d,    only: elarrays_2d, make_lmask
  use fftw,      only: dft
  implicit none

  private iu, pi, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qtt(nx,ny,D,rL,fC,T1,T2,glm,clm,gtype)
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T1, T2
  double complex, intent(out), dimension(nx,ny) :: glm, clm
  !(optional)
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
! [internal]
  integer :: i, j, nn(2)
  double precision, allocatable :: lmask(:,:), lx(:,:), ly(:,:), els(:,:), li(:,:)
  double complex, allocatable :: aT(:,:), almx(:,:), almy(:,:), blmx(:,:), blmy(:,:)

  nn = (/nx,ny/)

  allocate(lx(nx,ny),ly(nx,ny),li(nx,ny),els(nx,ny))
  call elarrays_2d(nn,D,elx=lx,ely=ly,eli=li,els=els)

  ! kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  ! filtering
  allocate(lmask(nx,ny))
  call make_lmask(nn,D,rL,lmask)

  allocate(aT(nx,ny),almx(nx,ny),almy(nx,ny),blmx(nx,ny),blmy(nx,ny))
  aT   = lmask*T1
  almx = 0.5d0*lmask*lx*fC*T2
  almy = 0.5d0*lmask*ly*fC*T2

  ! convolution
  call dft(aT,nn,D,-1)
  call dft(almx,nn,D,-1)
  call dft(almy,nn,D,-1)
  almx = aT*almx
  almy = aT*almy

  ! filtering
  aT   = lmask*T2
  blmx = 0.5d0*lmask*lx*fC*T1
  blmy = 0.5d0*lmask*ly*fC*T1

  ! convolution
  call dft(aT,nn,D,-1)
  call dft(blmx,nn,D,-1)
  call dft(blmy,nn,D,-1)
  blmx = aT*blmx
  blmy = aT*blmy

  deallocate(aT)

  ! to Fourier mode
  almx = almx + blmx
  almy = almy + blmy
  call dft(almx,nn,D,1)
  call dft(almy,nn,D,1)

  ! estimator 
  glm = (almx*lx+almy*ly)*li
  clm = (almx*ly-almy*lx)*li

  deallocate(almx,almy,blmx,blmy,els,lx,ly,li,lmask)

end subroutine qtt


subroutine qte(nx,ny,D,rL,fC,T,E,glm,clm,gtype)
  implicit none
! I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T, E
  double complex, intent(out), dimension(nx,ny) :: glm, clm
  !(optional)
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
! [internal]
  integer :: i, j, nn(2)
  double precision, allocatable :: lmask(:,:), lx(:,:), ly(:,:), els(:,:), li(:,:)
  double complex, allocatable :: aTx(:,:), aTy(:,:), bT(:,:), aE(:,:), bEx(:,:), bEy(:,:), almx(:,:), almy(:,:), ei2p(:,:)

  nn = (/nx,ny/)

  allocate(lx(nx,ny),ly(nx,ny),li(nx,ny),els(nx,ny),ei2p(nx,ny))
  call elarrays_2d(nn,D,elx=lx,ely=ly,eli=li,els=els,ei2p=ei2p)

  ! kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  ! filtering
  allocate(lmask(nx,ny))
  call make_lmask(nn,D,rL,lmask)

  allocate(aTx(nx,ny),aTy(nx,ny),bT(nx,ny),aE(nx,ny),bEx(nx,ny),bEy(nx,ny))

  aTx = lmask*lx*fC*T*ei2p
  aTy = lmask*ly*fC*T*ei2p
  aE  = lmask*E*conjg(ei2p)
  bT  = lmask*T
  bEx = lmask*lx*fC*E
  bEy = lmask*ly*fC*E

  deallocate(lmask,els,ei2p)

  ! convolution
  call dft(aE,nn,D,-1)
  call dft(bT,nn,D,-1)
  call dft(aTx,nn,D,-1)
  call dft(aTy,nn,D,-1)
  call dft(bEx,nn,D,-1)
  call dft(bEy,nn,D,-1)

  allocate(almx(nx,ny),almy(nx,ny))

  almx = iu*aimag(aE*aTx) + bT*bEx
  almy = iu*aimag(aE*aTy) + bT*bEy

  deallocate(aE,bEx,bEy,aTx,aTy,bT)

  call dft(almx,nn,D,1)
  call dft(almy,nn,D,1)

  ! form estimator 
  glm = (almx*lx+almy*ly)*li
  clm = (almx*ly-almy*lx)*li

  deallocate(almx,almy,lx,ly,li)

end subroutine qte


subroutine qtb(nx,ny,D,rL,fC,T,B,glm,clm,gtype)
  implicit none
! [input]
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: T, B
  double complex, intent(out), dimension(nx,ny) :: glm, clm
  !(optional)
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
! [internal]
  integer :: i, n, nn(2)
  double precision, allocatable :: lmask(:,:), lx(:,:), ly(:,:), els(:,:), li(:,:)
  double complex, allocatable :: wB(:,:),wTx(:,:),wTy(:,:),almx(:,:),almy(:,:),ei2p(:,:)

  nn   = (/nx,ny/)

  allocate(wB(nx,ny),wTx(nx,ny),wTy(nx,ny),lx(nx,ny),ly(nx,ny),li(nx,ny),els(nx,ny),ei2p(nx,ny))
  call elarrays_2d(nn,D,elx=lx,ely=ly,eli=li,els=els,ei2p=ei2p)

  ! kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  ! filtering
  allocate(lmask(nx,ny))
  call make_lmask(nn,D,rL,lmask)

  wB  = lmask*B*conjg(ei2p)
  wTx = lmask*lx*fC*T*ei2p
  wTy = lmask*ly*fC*T*ei2p

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wTx,nn,D,-1)
  call dft(wTy,nn,D,-1)

  allocate(almx(nx,ny),almy(nx,ny))
  almx = -iu*dble(wB*wTx)
  almy = -iu*dble(wB*wTy)

  deallocate(wB,wTx,wTy)
  call dft(almx,nn,D,1)
  call dft(almy,nn,D,1)

  ! form estimator 
  glm = (almx*lx+almy*ly)*li
  clm = (almx*ly-almy*lx)*li
  deallocate(almx,almy,lx,ly,li,els)

end subroutine qtb


subroutine qee(nx,ny,D,rL,fC,E1,E2,glm,clm,gtype)
  implicit none
! [input]
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: E1, E2
  double complex, intent(out), dimension(nx,ny) :: glm, clm
  !(optional)
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
! [internal]
  integer :: i, n, nn(2)
  double precision, allocatable :: lmask(:,:), lx(:,:), ly(:,:), els(:,:), li(:,:)
  double complex, allocatable :: wE1(:,:), wE2x(:,:), wE2y(:,:), almx(:,:), almy(:,:), ei2p(:,:)

  nn   = (/nx,ny/)

  allocate(wE1(nx,ny),wE2x(nx,ny),wE2y(nx,ny),lx(nx,ny),ly(nx,ny),li(nx,ny),els(nx,ny),ei2p(nx,ny))
  call elarrays_2d(nn,D,elx=lx,ely=ly,eli=li,els=els,ei2p=ei2p)

  ! kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  ! filtering
  allocate(lmask(nx,ny))
  call make_lmask(nn,D,rL,lmask)

  wE1  = lmask*E1*conjg(ei2p)
  wE2x = lmask*lx*fC*E2*ei2p
  wE2y = lmask*ly*fC*E2*ei2p

  ! convolution
  call dft(wE1,nn,D,-1)
  call dft(wE2x,nn,D,-1)
  call dft(wE2y,nn,D,-1)

  allocate(almx(nx,ny),almy(nx,ny))
  almx = iu*aimag(wE1*wE2x)
  almy = iu*aimag(wE1*wE2y)
  deallocate(wE1,wE2x,wE2y)

  call dft(almx,nn,D,1)
  call dft(almy,nn,D,1)

  ! estimator 
  glm = (lx*almx+ly*almy)*li
  clm = (ly*almy-lx*almy)*li
  deallocate(almx,almy,lx,ly,els,li)

end subroutine qee


subroutine qeb(nx,ny,D,rL,fC,E,B,glm,clm,gtype)
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: E, B
  double complex, intent(out), dimension(nx,ny) :: glm, clm
  !(optional)
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: i, n, nn(2)
  double precision, allocatable :: lmask(:,:), lx(:,:), ly(:,:), els(:,:), li(:,:)
  double complex, allocatable :: wB(:,:), almx(:,:), almy(:,:), wEx(:,:), wEy(:,:), ei2p(:,:)

  nn   = (/nx,ny/)

  allocate(wB(nx,ny),wEx(nx,ny),wEy(nx,ny),ei2p(nx,ny),lx(nx,ny),ly(nx,ny),els(nx,ny),li(nx,ny))
  call elarrays_2d(nn,D,elx=lx,ely=ly,els=els,eli=li,ei2p=ei2p)

  ! kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  ! filtering
  allocate(lmask(nx,ny))
  call make_lmask(nn,D,rL,lmask)

  wB  = lmask*B*conjg(ei2p)
  wEx = lmask*lx*fC*E*ei2p
  wEy = lmask*ly*fC*E*ei2p

  deallocate(ei2p)

  ! convolution
  call dft(wB,nn,D,-1)
  call dft(wEx,nn,D,-1)
  call dft(wEy,nn,D,-1)

  allocate(almx(nx,ny),almy(nx,ny))
  almx = -iu*dble(wB*wEx)
  almy = -iu*dble(wB*wEy)
  deallocate(wB,wEx,wEy)

  call dft(almx,nn,D,1)
  call dft(almy,nn,D,1)

  ! estimator 
  glm = (lx*almx+ly*almy) * li
  clm = (ly*almx-lx*almy) * li

  deallocate(almx,almy,lx,ly,els,li)

end subroutine qeb


subroutine qbb(nx,ny,D,rL,fC,B1,B2,glm,clm,gtype)
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: fC
  double complex, intent(in), dimension(nx,ny) :: B1, B2
  double complex, intent(out), dimension(nx,ny) :: glm, clm
  !(optional)
  character(1), intent(in), optional :: gtype
  !f2py character(1) :: gtype = ''
  !internal
  integer :: i, n, nn(2)
  double precision, allocatable :: lmask(:,:), lx(:,:), ly(:,:), els(:,:), li(:,:)
  double complex, allocatable :: wB1(:,:), wB2x(:,:), wB2y(:,:), almx(:,:), almy(:,:), ei2p(:,:)

  nn   = (/nx,ny/)

  allocate(wB1(nx,ny),wB2x(nx,ny),wB2y(nx,ny),ei2p(nx,ny),lx(nx,ny),ly(nx,ny),els(nx,ny),li(nx,ny))
  call elarrays_2d(nn,D,elx=lx,ely=ly,els=els,ei2p=ei2p)

  ! kappa factor
  li = -2d0*li**2
  if(.not.present(gtype).or.gtype/='k') li = 1d0

  ! filtering
  allocate(lmask(nx,ny))
  call make_lmask(nn,D,rL,lmask)

  allocate(wB1(nx,ny),wB2x(nx,ny),wB2y(nx,ny))
  wB1  = lmask*B1*conjg(ei2p)
  wB2x = lmask*lx*fC*B2*ei2p
  wB2y = lmask*ly*fC*B2*ei2p

  ! convolution
  call dft(wB1,nn,D,-1)
  call dft(wB2x,nn,D,-1)
  call dft(wB2y,nn,D,-1)

  allocate(almx(nx,ny),almy(nx,ny))
  almx = iu*aimag(wB1*wB2x)
  almy = iu*aimag(wB1*wB2y)

  deallocate(wB1,wB2x,wB2y)
  call dft(almx,nn,D,1)
  call dft(almy,nn,D,1)

  ! form estimator 
  glm = (almx*lx+almy*ly)*li
  clm = (almx*ly-almy*lx)*li
  deallocate(almx,almy)

end subroutine qbb


end module rec_lens


