!////////////////////////////////////////////////////!
! Point source reconstruction in Fullsky
!////////////////////////////////////////////////////!

module rec_src
  use alm_tools, only: alm2map, map2alm
  use constants, only: iu

  private alm2map, map2alm
  private iu

contains 


subroutine qtt(lmax,rlmin,rlmax,Tlm1,Tlm2,slm,nside_t,verbose)
!*  Reconstructing point sources from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of output point-source alms
!*    :rlmin/rlmax (int)  : Minimum/Maximum multipole of CMB for reconstruction
!*    :Tlm1 [l,m] (dcmplx): 1st inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*    :Tlm2 [l,m] (dcmplx): 2nd inverse-variance filtered temperature alm, with bounds (0:rlmax,0:rlmax)
!*
!*  Args(optional):
!*    :nside_t (int)      : Nside for the convolution calculation
!*    :verbose (bool)     : Output messages, default to False
!*
!*  Returns:
!*    :slm [l,m] (dcmplx) : Point-source alm, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm1, Tlm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: slm
  !internal
  integer :: l, npix, nside
  double precision, allocatable :: map(:,:)
  double complex, allocatable :: alm(:,:,:)
  !opt4py :: nside_t = 0
  !opt4py :: verbose = False

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc src-TT estimator, nside = ', nside
  npix = 12*nside**2

  ! alm to map 
  allocate(alm(2,0:rlmax,0:rlmax)); alm = 0d0
  do l = rlmin, rlmax
    alm(1,l,0:l) = Tlm1(l,0:l)
    alm(2,l,0:l) = Tlm2(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rlmax,rlmax,alm(1:1,:,:),map(1,:))
  call alm2map(nside,rlmax,rlmax,alm(2:2,:,:),map(2,:))
  deallocate(alm)

  ! map to alm
  allocate(alm(1,0:lmax,0:lmax))
  call map2alm(nside,lmax,lmax,map(1,:)*map(2,:),alm)
  slm = 0d0
  do l = 1, lmax
    slm(l,0:l) = 0.5d0*alm(1,l,0:l)
  end do
  deallocate(map,alm)

end subroutine qtt


end module rec_src


