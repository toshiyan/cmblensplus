!////////////////////////////////////////////////////!
! * CMB pol rotation angle reconstruction in curvedsky
!////////////////////////////////////////////////////!

module rec_rot
  use alm_tools, only: alm2map, alm2map_spin, map2alm
  use constants, only: iu

  private alm2map, alm2map_spin, map2alm
  private iu

contains 


subroutine qeb(lmax,rlmin,rlmax,fCE,Elm,Blm,alm,nside,verbose)
!*  Reconstructing pol rotation angle from the EB quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output lensing potential alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :fCE [l] (double)   : EE spectrum, with bounds (0:rlmax)
!*    :Elm [l,m] (dcmplx) : Inverse-variance filtered E-mode alm, with bounds (0:rlmax,0:rlmax)
!*    :Blm [l,m] (dcmplx) : Inverse-variance filtered B-mode alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside (int)        : Nside for the convolution calculation, default to lmax
!*    :verbose (bool)     : Output messages, default to False
!*
!*  Returns:
!*    :alm [l,m] (dcmplx) : Rotation angle alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside
  double precision, intent(in), dimension(0:rlmax) :: fCE
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Elm, Blm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: alm
  !internal
  integer :: l, ns, npix
  double precision, allocatable :: map(:), A(:,:), A2(:,:)
  double complex, allocatable :: xlm(:,:,:)
  !opt4py :: nside = 0
  !opt4py :: verbose = False

  ns = nside
  if (nside==0)  ns = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)   write(*,*) 'calc rot-EB estimator with nside=', ns
  npix = 12*ns**2

  ! convolution
  allocate(xlm(2,0:rlmax,0:rlmax),A(0:npix-1,2),A2(0:npix-1,2))
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(2,l,0:l) = Blm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A)
  xlm = 0d0
  do l = rlmin, rlmax
    xlm(1,l,0:l) = fCE(l)*Elm(l,0:l)
  end do 
  call alm2map_spin(nside,rlmax,rlmax,2,xlm,A2)
  deallocate(xlm)

  allocate(map(0:npix-1))
  map = A(:,1)*A2(:,2)-A(:,2)*A2(:,1)
  deallocate(A,A2)

  allocate(xlm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map,xlm)
  deallocate(map)

  ! compute alm
  write(*,*) 'compute polarization rotation'
  alm = 0d0
  do l = 1, lmax
    alm(l,0:l) = -2d0*xlm(1,l,0:l)
  end do
  deallocate(xlm)

end subroutine qeb


end module rec_rot


