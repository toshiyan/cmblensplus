!///////////////////////////////////////////////////////////!
! 2D Lensing Reconstruction Normalization
!///////////////////////////////////////////////////////////!

module norm_lens
  use constants, only: iu, i1, pi, dlc, twopi
  use grid2d, only: elarrays_2d, make_lmask
  use fftw,  only: dft
  implicit none

  private iu, i1, pi, dlc, twopi
  private elarrays_2d, make_lmask
  private dft

contains 


subroutine qtt(nx,ny,D,rL,OT,TT,eL,Ag,Ac)
!*  Normalization of the temperature quadratic estimator for CMB lensing potential and its curl mode
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OT[lx,ly] (double) : Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :TT[lx,ly] (double) : Theory temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :Ag[lx,ly] (dcmplx) : Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
!*    :Ac[lx,ly] (dcmplx) : Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, eL, rL, D, OT, TT
  !f2py intent(out) Ag, Ac
  !f2py depend(nx) OT, TT, Ag, Ac
  !f2py depend(ny) OT, TT, Ag, Ac
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: OT, TT
  double precision, intent(out), dimension(nx,ny) :: Ag, Ac
  !internal
  integer :: i, j, nn(2)
  double precision, allocatable, dimension(:,:) :: lmask, lx, ly, Agg, Acc
  double complex, allocatable, dimension(:,:) :: A, Axx, Axy, Ayy, Bx, By

  nn = (/nx,ny/)

  allocate(lx(nx,ny),ly(nx,ny),Agg(nx,ny),Acc(nx,ny),A(nx,ny),Axx(nx,ny),Axy(nx,ny),Ayy(nx,ny),Bx(nx,ny),By(nx,ny),lmask(nx,ny))

  call elarrays_2d(nn,D,elx=lx,ely=ly)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  
  Axx = lmask * lx**2 * TT**2*OT
  Axy = lmask * lx*ly * TT**2*OT
  Ayy = lmask * ly**2 * TT**2*OT
  A   = lmask * OT
  Bx  = lmask * lx * TT * OT
  By  = lmask * ly * TT * OT

  ! convolution
  call dft(Axx,nn,D,-1)
  call dft(Axy,nn,D,-1)
  call dft(Ayy,nn,D,-1)
  call dft(A,nn,D,-1)
  call dft(Bx,nn,D,-1)
  call dft(By,nn,D,-1)
  Axx = A*Axx + Bx**2
  Axy = A*Axy + Bx*By
  Ayy = A*Ayy + By**2

  call dft(Axx,nn,D,1)
  call dft(Axy,nn,D,1)
  call dft(Ayy,nn,D,1)

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Agg = lmask * (lx**2*Axx + 2*lx*ly*Axy + ly**2*Ayy)
  Acc = lmask * (ly**2*Axx - 2*lx*ly*Axy + lx**2*Ayy)

  ! inversion
  Ag = 0d0
  Ac = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Agg(i,j)>0)  Ag(i,j) = 1d0/Agg(i,j)
      if (Acc(i,j)>0)  Ac(i,j) = 1d0/Acc(i,j)
    end do
  end do
  
  deallocate(lx,ly,Agg,Acc,A,Axx,Axy,Ayy,Bx,By,lmask)

end subroutine qtt

subroutine qte(nx,ny,D,rL,OT,OE,TE,eL,Ag,Ac)
!*  Normalization of the TE quadratic estimator for CMB lensing potential and its curl mode
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OT[lx,ly] (double) : Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :OE[lx,ly] (double) : Inverse of Observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :TE[lx,ly] (double) : Theory TE cross spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :Ag[lx,ly] (dcmplx) : Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
!*    :Ac[lx,ly] (dcmplx) : Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, eL, rL, D, TE, OT, OE
  !f2py intent(out) Ag, Ac
  !f2py depend(nx) TE, OT, OE, Ag, Ac
  !f2py depend(ny) TE, OT, OE, Ag, Ac
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: TE, OT, OE
  double precision, intent(out), dimension(nx,ny) :: Ag, Ac
  !internal
  integer :: i, j, nn(2)
  double precision, allocatable, dimension(:,:) :: lmask, lx, ly, Agg, Acc
  double complex, allocatable, dimension(:,:) :: ei2p, Axx, Axy, Ayy, Cxx, Cxy, Cyy, Ax, Ay, A, Bm, Bp, Bx, By, Bxx, Bxy, Byy

  nn = (/nx,ny/)
  
  allocate(lmask(nx,ny),lx(nx,ny),ly(nx,ny),Agg(nx,ny),Acc(nx,ny))
  allocate(ei2p(nx,ny),Axx(nx,ny),Axy(nx,ny),Ayy(nx,ny),Cxx(nx,ny),Cxy(nx,ny),Cyy(nx,ny),Ax(nx,ny),Ay(nx,ny),A(nx,ny),Bm(nx,ny),Bp(nx,ny),Bx(nx,ny),By(nx,ny),Bxx(nx,ny),Bxy(nx,ny),Byy(nx,ny))

  call elarrays_2d(nn,D,lx,ly,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)

  Axx = lmask * OT * TE**2 * ei2p**2/2d0 * lx**2
  Axy = lmask * OT * TE**2 * ei2p**2/2d0 * lx*ly*2d0
  Ayy = lmask * OT * TE**2 * ei2p**2/2d0 * ly**2
  Bm  = lmask * OE * conjg(ei2p)**2
  Cxx = lmask * OT * TE**2 / 2d0 * lx**2
  Cxy = lmask * OT * TE**2 / 2d0 * lx*ly*2d0
  Cyy = lmask * OT * TE**2 / 2d0 * ly**2
  Bp  = lmask * OE
  Ax  = lmask * OT * TE * ei2p * lx
  Ay  = lmask * OT * TE * ei2p * ly
  Bx  = lmask * OE * TE * conjg(ei2p) * lx
  By  = lmask * OE * TE * conjg(ei2p) * ly
  A   = lmask * OT
  Bxx = lmask * OE * TE**2 * lx**2
  Bxy = lmask * OE * TE**2 * lx*ly*2d0
  Byy = lmask * OE * TE**2 * ly**2

  ! inverse FT
  call dft(Axx,nn,D,-1)
  call dft(Axy,nn,D,-1)
  call dft(Ayy,nn,D,-1)
  call dft(Cxx,nn,D,-1)
  call dft(Cxy,nn,D,-1)
  call dft(Cyy,nn,D,-1)
  call dft(Ax,nn,D,-1)
  call dft(Ay,nn,D,-1)
  call dft(A,nn,D,-1)
  call dft(Bm,nn,D,-1)
  call dft(Bp,nn,D,-1)
  call dft(Bx,nn,D,-1)
  call dft(By,nn,D,-1)
  call dft(Bxx,nn,D,-1)
  call dft(Bxy,nn,D,-1)
  call dft(Byy,nn,D,-1)

  ! convolution
  Axx = Axx*Bm + Cxx*Bp + 2*Ax*Bx + A*Bxx
  Axy = Axy*Bm + Cxy*Bp + 2*(Ax*By + Ay*Bx) + A*Bxy
  Ayy = Ayy*Bm + Cyy*Bp + 2*Ay*By + A*Byy
  call dft(Axx,nn,D,1)
  call dft(Axy,nn,D,1)
  call dft(Ayy,nn,D,1)

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Agg = lmask * (lx**2*Axx + lx*ly*Axy + ly**2*Ayy)
  Acc = lmask * (ly**2*Axx - lx*ly*Axy + lx**2*Ayy)

  ! inversion
  Ag = 0d0
  Ac = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Agg(i,j)>0)  Ag(i,j) = 1d0/Agg(i,j)
      if (Acc(i,j)>0)  Ac(i,j) = 1d0/Acc(i,j)
    end do
  end do
  
  deallocate(lmask, lx, ly, Agg, Acc, ei2p, Axx, Axy, Ayy, Cxx, Cxy, Cyy, Ax, Ay, A, Bm, Bp, Bx, By, Bxx, Bxy, Byy)


end subroutine qte

subroutine qtb(nx,ny,D,OT,OB,TE,rL,eL,Ag,Ac)
!*  Normalization of the TB quadratic estimator for CMB lensing potential and its curl mode
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OT[lx,ly] (double) : Inverse of Observed temperature power spectrum on 2D grid, with bounds (nx,ny)
!*    :OB[lx,ly] (double) : Inverse of Observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :TE[lx,ly] (double) : Theory TE cross spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :Ag[lx,ly] (dcmplx) : Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
!*    :Ac[lx,ly] (dcmplx) : Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, eL, rL, D, TE, OT, OB
  !f2py intent(out) Ag, Ac
  !f2py depend(nx) TE, OT, OB, Ag, Ac
  !f2py depend(ny) TE, OT, OB, Ag, Ac
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: TE, OT, OB
  double precision, intent(out), dimension(nx,ny) :: Ag, Ac
  !internal
  integer :: i, j, nn(2)
  double precision, allocatable, dimension(:,:) :: lx, ly, lmask, Agg, Acc
  double complex, allocatable, dimension(:,:) :: Axx, Axy, Ayy, Bxx, Bxy, Byy, A, B, ei2p

  nn = (/nx,ny/)
  
  allocate(lx(nx,ny),ly(nx,ny),lmask(nx,ny),Agg(nx,ny),Acc(nx,ny),Axx(nx,ny),Axy(nx,ny),Ayy(nx,ny),Bxx(nx,ny),Bxy(nx,ny),Byy(nx,ny),A(nx,ny),B(nx,ny),ei2p(nx,ny))
  
  call elarrays_2d(nn,D,lx,ly,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  Bxx = lmask * lx**2*TE**2*OT
  Bxy = lmask * 2*lx*ly*TE**2*OT
  Byy = lmask * ly**2*TE**2*OT
  B   = lmask * 0.5d0*OB
  Axx = -ei2p**2 * Bxx
  Axy = -ei2p**2 * Bxy
  Ayy = -ei2p**2 * Byy
  A   = conjg(ei2p)**2 * B

  ! convolution
  call dft(A,nn,D,-1)
  call dft(B,nn,D,-1)
  call dft(Axx,nn,D,-1)
  call dft(Axy,nn,D,-1)
  call dft(Ayy,nn,D,-1)
  call dft(Bxx,nn,D,-1)
  call dft(Bxy,nn,D,-1)
  call dft(Byy,nn,D,-1)
  Axx = A*Axx + B*Bxx
  Axy = A*Axy + B*Bxy
  Ayy = A*Ayy + B*Byy
  call dft(Axx,nn,D,1)
  call dft(Axy,nn,D,1)
  call dft(Ayy,nn,D,1)
 
  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Agg = lmask * (lx**2*Axx + lx*ly*Axy + ly**2*Ayy)
  Acc = lmask * (ly**2*Axx - lx*ly*Axy + lx**2*Ayy)

  ! inversion
  Ag = 0d0
  Ac = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Agg(i,j)>0)  Ag(i,j) = 1d0/Agg(i,j)
      if (Acc(i,j)>0)  Ac(i,j) = 1d0/Acc(i,j)
    end do
  end do

  deallocate(lx,ly,lmask,Agg,Acc,Axx,Axy,Ayy,Bxx,Bxy,Byy,A,B,ei2p)

end subroutine qtb

subroutine qee(nx,ny,D,OE,EE,rL,eL,Ag,Ac)
!*  Normalization of the EE quadratic estimator for CMB lensing potential and its curl mode
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OE[lx,ly] (double) : Inverse of Observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :EE[lx,ly] (double) : Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :Ag[lx,ly] (dcmplx) : Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
!*    :Ac[lx,ly] (dcmplx) : Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, eL, rL, D, EE, OE
  !f2py intent(out) Ag, Ac
  !f2py depend(nx) EE, OE, Ag, Ac
  !f2py depend(ny) EE, OE, Ag, Ac
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: EE, OE
  double precision, intent(out), dimension(nx,ny) :: Ag, Ac
! [internal]
  integer :: i, j, nn(2)
  double precision, allocatable, dimension(:,:) :: lx, ly, lmask, Agg, Acc
  double complex, allocatable, dimension(:,:) :: Axx, Axy, Ayy, Bxx, Bxy, Byy, A, B, Ax, Ay, Bx, By, ei2p

  nn = (/nx,ny/)
  
  allocate(lx(nx,ny),ly(nx,ny),lmask(nx,ny),Agg(nx,ny),Acc(nx,ny),Axx(nx,ny),Axy(nx,ny),Ayy(nx,ny),Bxx(nx,ny),Bxy(nx,ny),Byy(nx,ny),A(nx,ny),B(nx,ny),Ax(nx,ny),Ay(nx,ny),Bx(nx,ny),By(nx,ny),ei2p(nx,ny))
  
  call elarrays_2d(nn,D,lx,ly,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)

  Bxx = lmask * (lx**2*EE**2*OE)
  Bxy = lmask * (2*lx*ly*EE**2*OE)
  Byy = lmask * (ly**2*EE**2*OE)
  B   = lmask * OE
  Axx = ei2p**2 * Bxx
  Axy = ei2p**2 * Bxy
  Ayy = ei2p**2 * Byy
  A   = B * conjg(ei2p)**2

  Bx  = lmask * lx*EE*OE
  By  = lmask * ly*EE*OE
  Ax  = Bx * ei2p**2
  Ay  = By * ei2p**2
  
  ! convolution
  call dft(A,nn,D,-1)
  call dft(B,nn,D,-1)
  call dft(Axx,nn,D,-1)
  call dft(Axy,nn,D,-1)
  call dft(Ayy,nn,D,-1)
  call dft(Bxx,nn,D,-1)
  call dft(Bxy,nn,D,-1)
  call dft(Byy,nn,D,-1)
  call dft(Ax,nn,D,-1)
  call dft(Ay,nn,D,-1)
  call dft(Bx,nn,D,-1)
  call dft(By,nn,D,-1)

  Axx = A*Axx - Ax*conjg(Ax) + B*Bxx + Bx**2
  Axy = A*Axy - Ax*conjg(Ay) - Ay*conjg(Ax) + B*Bxy + 2*Bx*By
  Ayy = A*Ayy - Ay*conjg(Ay) + B*Byy + By**2

  call dft(Axx,nn,D,1)
  call dft(Axy,nn,D,1)
  call dft(Ayy,nn,D,1)

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Agg = lmask * (lx**2*Axx + lx*ly*Axy + ly**2*Ayy)
  Acc = lmask * (ly**2*Axx - lx*ly*Axy + lx**2*Ayy)

  ! inversion
  Ag = 0d0
  Ac = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Agg(i,j)>0)  Ag(i,j) = 2d0/Agg(i,j)
      if (Acc(i,j)>0)  Ac(i,j) = 2d0/Acc(i,j)
    end do
  end do
  
  deallocate(lx,ly,lmask,Agg,Acc,Axx,Axy,Ayy,Bxx,Bxy,Byy,A,B,Ax,Ay,Bx,By,ei2p)

end subroutine qee

subroutine qeb(nx,ny,D,OE,OB,EE,rL,eL,Ag,Ac)
!*  Normalization of the EB quadratic estimator for CMB lensing potential and its curl mode
!*
!*  Args:
!*    :nx, ny (int)       : Number of Lx and Ly grids
!*    :D[xy] (double)     : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :rL[2] (int)        : Minimum and maximum multipole of CMB for reconstruction
!*    :OE[lx,ly] (double) : Inverse of Observed E-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :OB[lx,ly] (double) : Inverse of Observed B-mode power spectrum on 2D grid, with bounds (nx,ny)
!*    :EE[lx,ly] (double) : Theory E-mode spectrum on 2D grid, with bounds (nx,ny)
!*    :eL[2] (int)        : Minimum and maximum multipole of output normalization spectrum, with bounds (2)
!*
!*  Returns:
!*    :Ag[lx,ly] (dcmplx) : Normalization of CMB lensing potential on 2D grid, with bounds (nx,ny)
!*    :Ac[lx,ly] (dcmplx) : Normalization of Curl mode (pseudo lensing potential) on 2D grid, with bounds (nx,ny)
!*
  !f2py intent(in) nx, ny, eL, rL, D, EE, OE, OB
  !f2py intent(out) Ag, Ac
  !f2py depend(nx) EE, OE, OB, Ag, Ac
  !f2py depend(ny) EE, OE, OB, Ag, Ac
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: EE, OE, OB
  double precision, intent(out), dimension(nx,ny) :: Ag, Ac
  !internal
  integer :: i, j, nn(2)
  double precision, allocatable, dimension(:,:) :: lx, ly, lmask, Agg, Acc
  double complex, allocatable, dimension(:,:) :: Axx, Axy, Ayy, Bxx, Bxy, Byy, A, B, ei2p

  nn = (/nx,ny/)
  
  allocate(lx(nx,ny),ly(nx,ny),lmask(nx,ny),Agg(nx,ny),Acc(nx,ny),Axx(nx,ny),Axy(nx,ny),Ayy(nx,ny),Bxx(nx,ny),Bxy(nx,ny),Byy(nx,ny),A(nx,ny),B(nx,ny),ei2p(nx,ny))

  call elarrays_2d(nn,D,lx,ly,ei2p=ei2p)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  Bxx = lmask * lx**2*EE**2*OE
  Bxy = lmask * 2*lx*ly*EE**2*OE
  Byy = lmask * ly**2*EE**2*OE
  B   = lmask * 0.5d0*OB
  Axx = -ei2p**2 * Bxx
  Axy = -ei2p**2 * Bxy
  Ayy = -ei2p**2 * Byy
  A   = B * conjg(ei2p)**2

  ! convolution
  call dft(A,nn,D,-1)
  call dft(B,nn,D,-1)
  call dft(Axx,nn,D,-1)
  call dft(Axy,nn,D,-1)
  call dft(Ayy,nn,D,-1)
  call dft(Bxx,nn,D,-1)
  call dft(Bxy,nn,D,-1)
  call dft(Byy,nn,D,-1)
  Axx = A*Axx + B*Bxx
  Axy = A*Axy + B*Bxy
  Ayy = A*Ayy + B*Byy
  call dft(Axx,nn,D,1)
  call dft(Axy,nn,D,1)
  call dft(Ayy,nn,D,1)

  ! normalization
  call make_lmask(nn,D,eL,lmask)
  Agg = lmask * (lx**2*Axx + lx*ly*Axy + ly**2*Ayy)
  Acc = lmask * (ly**2*Axx - lx*ly*Axy + lx**2*Ayy)

  ! inversion
  Ag = 0d0
  Ac = 0d0
  do i = 1, nx
    do j = 1, ny
      if (Agg(i,j)>0)  Ag(i,j) = 1d0/Agg(i,j)
      if (Acc(i,j)>0)  Ac(i,j) = 1d0/Acc(i,j)
    end do
  end do

  deallocate(lx,ly,lmask,Agg,Acc,Axx,Axy,Ayy,Bxx,Bxy,Byy,A,B,ei2p)

end subroutine qeb



end module norm_lens



