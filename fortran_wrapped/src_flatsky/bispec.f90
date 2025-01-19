!////////////////////////////////////////////////////!
! bispectrum measurement from flatsky simulations
!////////////////////////////////////////////////////!

module bispec
  use constants, only: dlc
  use grid2d,    only: make_binmask
  use fftw,      only: dft
  use fft_utils, only: prep_filtered_map

  private dlc
  private make_binmask
  private dft
  private prep_filtered_map

contains

subroutine bispec_norm(nx,ny,D,bp,dbin_max,bn,norm)
  !f2py intent(in) nx, ny, dbin_max, bn, D, bp
  !f2py intent(out) norm
  !f2py depend(bn) bp, norm
  implicit none
  !I/O
  integer, intent(in) :: nx, ny, dbin_max, bn
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(bn+1) :: bp
  double precision, intent(out), dimension(bn,bn,bn) :: norm
  !internal
  double complex, allocatable :: ulm(:,:,:)
  double precision, allocatable :: kmap(:,:,:,:)
  !opt4py :: bn = 1
  !add2py :: bn = len(bp) - 1
  !opt4py :: dbin_max = -1
  !add2py :: if dbmax==-1: dbin_max = bn

  allocate(ulm(1,nx,ny),kmap(1,bn,nx,ny)); kmap=0d0

  ! first prepare filtered map for each bin
  ulm = 1d0
  call prep_filtered_map(nx,ny,D,bp,ulm,kmap)

  ! next compute the product of filtered maps for each combination of b1, b2, b3
  call bispec_bin(1,bn,nx,ny,kmap,bp,dbin_max,norm)

  norm = norm*D(1)*D(2)

  deallocate(ulm,kmap)

end subroutine bispec_norm

subroutine bispec_bin(kn,bn,nx,ny,kmap,bp,dbin_max,bispec)
  !unnormalized binned bispectrum for 3D shape
  !f2py intent(in) kn, bn, nx, ny, dbin_max, kmap!multipole-filteredmap, bp
  !f2py intent(out) bispec
  !f2py depend(kn) kmap!multipole-filteredmap
  !f2py depend(bn) kmap!multipole-filteredmap, bp, bispec
  !f2py depend(nx) kmap!multipole-filteredmap
  !f2py depend(ny) kmap!multipole-filteredmap
  implicit none
  !I/O
  integer, intent(in) :: kn, bn, nx, ny, dbin_max
  double precision, intent(in), dimension(kn,bn,nx,ny) :: kmap  !multipole-filtered map
  double precision, intent(in), dimension(bn+1) :: bp
  double precision, intent(out), dimension(bn,bn,bn) :: bispec
  !internal
  integer :: b0, b1, b2, bmax
  double precision :: lbmin, lbmax
  double precision :: Bb111, Bb112, Bb121, Bb211

  !opt4py :: bn = 1
  !add2py :: bn = len(bp) - 1
  !opt4py :: kn = 1
  !add2py :: kn = len(kmap[:,0,0,0])
  !opt4py :: nx = 0
  !add2py :: nx = len(kmap[0,0,:,0])
  !opt4py :: ny = 0
  !add2py :: ny = len(kmap[0,0,0,:])
  !opt4py :: dbin_max = -1
  !add2py :: if dbmax==-1: dbin_max = bn

  bispec = 0d0

  do b0 = 1, bn
    do b1 = b0, bn
      if (b1>b0+dbin_max) cycle !discard far separate bins 
      !assuming bp(b0)<bp(b1)
      lbmin = bp(b1) - bp(b0+1)
      lbmax = bp(b1+1) + bp(b0+1)
      do b2 = b1, bn
        if (b2>b1+dbin_max) cycle !discard far separate bins
        if (bp(b2)<lbmin.or.lbmax<bp(b2+1)) cycle !triangle condition
        select case(kn)
        case(1)
          Bb111 = sum(kmap(1,b0,:,:)*kmap(1,b1,:,:)*kmap(1,b2,:,:))
          bispec(b0,b1,b2) = Bb111
        case(2)
          Bb112 = sum(kmap(1,b0,:,:)*kmap(1,b1,:,:)*kmap(2,b2,:,:))
          Bb121 = sum(kmap(1,b0,:,:)*kmap(2,b1,:,:)*kmap(1,b2,:,:))
          Bb211 = sum(kmap(2,b0,:,:)*kmap(1,b1,:,:)*kmap(1,b2,:,:))
          bispec(b0,b1,b2) = Bb112 + Bb121 + Bb211
        case default
          stop 'only kn=1 or 2 types are supported'
        end select
      end do
    end do
  end do

end subroutine bispec_bin

subroutine binfilter(nx,ny,D,bp,bf,bn)
!*  The multipole bin binary-mask given on the 2D Fourier grids so that M = 1 inside the multipole bin and 0 otherwise
!*
!*  Args:
!*    :nx, ny (int)          : Number of Lx and Ly grids
!*    :D[2] (double)         : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :bp[bn+1] (double)     : Multipole bin edges
!*
!*  Returns:
!*    :bf[bn,nx,ny] (double) : The multipole bin binary-mask for each multipol bin

  !f2py intent(in) nx, ny, bn, D, bp
  !f2py intent(out) bf
  !f2py depend(bn) bp, bf
  !f2py depend(nx) bf
  !f2py depend(ny) bf
  implicit none
  integer, intent(in) :: nx, ny, bn
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(bn+1) :: bp
  double precision, intent(out), dimension(bn,nx,ny) :: bf
  !opt4py :: bn = 1
  !add2py :: bn = len(bp) - 1
  
  call make_binmask((/nx,ny/),D,bp(1:bn),bp(2:bn+1),bf,bn)

end subroutine binfilter

subroutine bispec_norm_1d(nx,ny,D,bfs,bnorm,bn)
!*  Normalization of the 1D binned bispectrum estimator
!*
!*  Args:
!*    :nx, ny (int)             : Number of Lx and Ly grids
!*    :D[2] (double)            : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :bfs[3,bn,nx,ny] (double) : Multipole bin mask on 2D grids obtained from the binfilter function
!*
!*  Returns:
!*    :bnorm[bn] (double)       : Normalization of 1D binned bispectrum at each multipole bin

  !f2py intent(in) nx, ny, bn, D, bfs
  !f2py intent(out) bnorm
  !f2py depend(nx) bfs
  !f2py depend(ny) bfs
  !f2py depend(bn) bfs, bnorm
  implicit none
  !I/O
  integer, intent(in) :: nx, ny, bn
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(3,bn,nx,ny) :: bfs
  double precision, intent(out), dimension(bn) :: bnorm
  !internal
  integer :: b, i, nn(2)
  complex(dlc), allocatable :: wlm(:,:,:)
  
  !opt4py :: bn = 1
  !add2py :: bn = len(bfs[0,:,0,0])

  nn = (/nx,ny/)
  bnorm = 0d0

  do b = 1, bn
    allocate(wlm(3,nx,ny));  wlm=0d0
    do i = 1, 3
      wlm(i,:,:) = bfs(i,b,:,:)
      call dft(wlm(i,:,:),nn,D,-1)
    end do
    bnorm(b) = D(1)*D(2)*sum(wlm(1,:,:)*wlm(2,:,:)*wlm(3,:,:))
    deallocate(wlm)
  end do

end subroutine bispec_norm_1d

subroutine bispec_bin_1d(nx,ny,D,bfs,bnorm,alm,bispec,bn)
!*  1D binned bispectrum estimator
!*
!*  Args:
!*    :nx, ny (int)             : Number of Lx and Ly grids
!*    :D[2] (double)            : Map side length, or equivalent to dLx/2pi, dLy/2pi, with bounds (2)
!*    :bfs[3,bn,nx,ny] (double) : Multipole bin mask on 2D grids obtained from the binfilter function
!*    :bnorm[bn] (double)       : Normalization of 1D binned bispectrum at each multipole bin
!*    :alm[3,nx,ny] (dcmplx)    : Fourier modes for each leg
!*
!*  Returns:
!*    :bispec[bn] (double)      : 1D binned bispectrum at each multipole bin

  !f2py intent(in) nx, ny, bn, D, bfs, bnorm, alm
  !f2py intent(out) bispec
  !f2py depend(nx) bfs, alm
  !f2py depend(ny) bfs, alm
  !f2py depend(bn) bfs, bnorm, bispec
  implicit none
  !I/O
  integer, intent(in) :: nx, ny, bn
  double precision, intent(in), dimension(2) :: D
  double precision, intent(in), dimension(3,bn,nx,ny) :: bfs
  double precision, intent(in), dimension(bn) :: bnorm
  double precision, intent(out), dimension(bn) :: bispec
  double complex, intent(in), dimension(3,nx,ny) :: alm
  integer :: b, i, nn(2)
  double complex, allocatable :: walm(:,:,:)

  !opt4py :: bn = 1
  !add2py :: bn = len(bfs[0,:,0,0])

  nn = (/nx,ny/)
  bispec = 0d0
  
  do b = 1, bn
    if (bnorm(b) == 0d0) cycle
    allocate(walm(3,nx,ny));  walm=0d0
    do i = 1, 3
      walm(i,:,:) = bfs(i,b,:,:)*alm(i,:,:)
      call dft(walm(i,:,:),nn,D,-1)
    end do
    bispec(b) = sum(walm(1,:,:)*walm(2,:,:)*walm(3,:,:))/bnorm(b)
    deallocate(walm)
  end do

end subroutine bispec_bin_1d



end module bispec


