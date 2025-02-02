!////////////////////////////////////////////////////!
! Amplitude modulation reconstruction
!////////////////////////////////////////////////////!

module rec_iamp
  use alm_tools, only: alm2map, map2alm, alm2map_spin
  use constants, only: iu

  private alm2map, map2alm, alm2map_spin
  private iu

contains 


subroutine qeb(lmax,rlmin,rlmax,EB,Elm,Blm,alm,nside_t,verbose)
!*  Reconstructing amplitude modulation by the odd EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :EB [l] (double)    : EB spectrum, with bounds (0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)      : Nside for the convolution calculation
!*    :verbose (bool)     : Output messages, default to False
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Reconstructed alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  double precision, intent(in), dimension(0:rlmax) :: EB
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, nside, npix
  double precision, allocatable :: map(:), A(:,:), A2(:,:)
  double complex, allocatable :: xlm(:,:,:)
  !opt4py :: nside_t = 0
  !opt4py :: verbose = False

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc tau-EB odd estimator with nside=', nside
  npix = 12*nside**2

  ! convolution 1
  allocate(xlm(2,0:rlmax,0:rlmax),A(0:npix-1,2),A2(0:npix-1,2))
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(1,l,0:l) = Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A)
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(2,l,0:l) = EB(l)*Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A2)

  deallocate(xlm)

  allocate(map(0:npix-1))
  map = -A(:,1)*A2(:,2)+A(:,2)*A2(:,1)
  deallocate(A,A2)

  ! convolution 2
  allocate(xlm(2,0:rlmax,0:rlmax),A(0:npix-1,2),A2(0:npix-1,2))
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(1,l,0:l) = EB(l)*Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A2)
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(2,l,0:l) = Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A)
  deallocate(xlm)

  map = map + A(:,1)*A2(:,2)-A(:,2)*A2(:,1)
  deallocate(A,A2)

  allocate(xlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map,xlm)
  deallocate(map)

  ! compute alm 
  write(*,*) 'compute reconstructed fields'
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = -2d0*xlm(1,l,0:l)
  end do
  deallocate(xlm)

end subroutine qeb


end module rec_iamp


