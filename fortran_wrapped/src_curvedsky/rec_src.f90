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
  !f2py intent(in) verbose, lmax, rlmin, rlmax, nside_t, Tlm1, Tlm2
  !f2py intent(out) slm
  !f2py depend(rlmax) Tlm1, Tlm2
  !f2py depend(lmax) slm
  !I/O
  logical, intent(in) :: verbose
  integer, intent(in) :: lmax, rlmin, rlmax, nside_t
  double complex, intent(in), dimension(0:rlmax,0:rlmax) :: Tlm1, Tlm2
  double complex, intent(out), dimension(0:lmax,0:lmax) :: slm
  !internal
  integer :: l, npix, nside
  double precision, allocatable :: map1(:), map2(:)
  double complex, allocatable :: alm1(:,:,:), alm2(:,:,:)
  !opt4py :: nside_t = 0
  !opt4py :: verbose = False

  nside = nside_t
  if (nside_t==0)  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))
  if (verbose)  write(*,*) 'calc src-TT estimator, nside = ', nside
  npix = 12*nside**2

  ! alm to map 
  allocate(alm1(1,0:rlmax,0:rlmax),alm2(1,0:rlmax,0:rlmax)); alm1=0d0; alm2=0d0
  do l = rlmin, rlmax
    alm1(1,l,0:l) = Tlm1(l,0:l)
    alm2(1,l,0:l) = Tlm2(l,0:l)
  end do 
  allocate(map1(0:npix-1),map2(0:npix-1))
  call alm2map(nside,rlmax,rlmax,alm1,map1)
  call alm2map(nside,rlmax,rlmax,alm2,map2)
  deallocate(alm1,alm2)

  ! map to alm
  allocate(alm1(1,0:lmax,0:lmax))
  map1 = map1*map2
  call map2alm(nside,lmax,lmax,map1,alm1)
  slm = 0d0
  do l = 1, lmax
    slm(l,0:l) = 0.5d0*alm1(1,l,0:l)
  end do
  deallocate(map1,map2,alm1)

end subroutine qtt



end module rec_src



