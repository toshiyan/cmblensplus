!///////////////////////////////////////////////////////////!
! 2D Lensing Reconstruction Normalization
!///////////////////////////////////////////////////////////!

module norm_lens
  use constants, only: iu, i1, pi, dlc, twopi
  use grid2d, only: elarrays_2d, make_lmask
  use fftw, only: dft
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
  implicit none
  !I/O
  integer, intent(in) :: nx, ny
  integer, intent(in), dimension(2) :: eL, rL
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(nx,ny) :: OT, TT
  double precision, intent(out), dimension(nx,ny) :: Ag, Ac
  !internal
  integer :: i, j, nn(2)
  double precision, allocatable, dimension(:,:) :: lmask, lx, ly, Agg, Acc, iAgg, iAcc
  complex(dlc), allocatable, dimension(:,:) :: A, Axx, Axy, Ayy, Bx, By

  nn = (/nx,ny/)

  allocate(lx(nx,ny),ly(nx,ny),Agg(nx,ny),Acc(nx,ny),iAgg(nx,ny),iAcc(nx,ny),A(nx,ny),Axx(nx,ny),Axy(nx,ny),Ayy(nx,ny),Bx(nx,ny),By(nx,ny),lmask(nx,ny))

  call elarrays_2d(nn,D,elx=lx,ely=ly)

  ! filtering
  call make_lmask(nn,D,rL,lmask)
  !lmask = 1d0
  write(*,*) 'X'
  
  Axx = lmask * lx**2 * TT**2*OT
  Axy = lmask * lx*ly * TT**2*OT
  Ayy = lmask * ly**2 * TT**2*OT
  A   = lmask * OT
  Bx  = lmask * lx * TT * OT
  By  = lmask * ly * TT * OT

  ! convolution
  write(*,*) 'G'
  call dft(Axx,nn,D,-1)
  write(*,*) 'F'
  !call dft(Axy,nn,D,-1)
  write(*,*) 'F'
  !call dft(Ayy,nn,D,-1)
  write(*,*) 'F'
  !call dft(A,nn,D,-1)
  write(*,*) 'F'
  !call dft(Bx,nn,D,-1)
  write(*,*) 'F'
  !call dft(By,nn,D,-1)
  write(*,*) 'F'
  Axx = A*Axx + Bx**2
  write(*,*) 'F'
  Axy = A*Axy + Bx*By
  write(*,*) 'F'
  Ayy = A*Ayy + By**2

  write(*,*) 'X'
  !call dft(Axx,nn,D,1)
  write(*,*) 'F'
  !call dft(Axy,nn,D,1)
  write(*,*) 'F'
  !call dft(Ayy,nn,D,1)

  write(*,*) 'F'
  ! normalization
  call make_lmask(nn,D,eL,lmask)
  write(*,*) 'B'
  iAgg = lmask * (lx**2*Axx + 2*lx*ly*Axy + ly**2*Ayy)
  iAcc = lmask * (ly**2*Axx - 2*lx*ly*Axy + lx**2*Ayy)

  ! inversion
  Agg = 0d0
  Acc = 0d0
  write(*,*) 'D'
  do i = 1, nx
    do j = 1, ny
      if (iAgg(i,j)>0)  Agg(i,j) = 1d0/iAgg(i,j)
      if (iAcc(i,j)>0)  Acc(i,j) = 1d0/iAcc(i,j)
    end do
  end do
  
  Ag = Agg
  Ac = Acc

  deallocate(lx,ly,Agg,Acc,iAgg,iAcc,A,Axx,Axy,Ayy,Bx,By,lmask)

end subroutine qtt



end module norm_lens


