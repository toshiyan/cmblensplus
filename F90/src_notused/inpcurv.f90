!//////////////////////////////////////////////////////////////////!

module inpcurv
  implicit none

contains 


subroutine mat_multi(npix,lmax,clh,nij,x,v,mtype)
!* multiplying matrix
  implicit none
  !I/O
  integer, intent(in) :: npix, lmax
  double precision, intent(in), dimension(0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(0:npix-1) :: nij
  double complex, intent(in), dimension(0:lmax,0:lmax) :: x
  double complex, intent(out), dimension(0:lmax,0:lmax) :: v
  !optional
  character(4), intent(in), optional :: mtype
  !internal
  character(4) :: m = ''
  integer :: l, nside
  double precision :: map(0:npix-1)
  double complex :: alm(1,0:lmax,0:lmax)

  if (present(mtype))  m = mtype

  nside = int(sqrt(npix/12d0))

  if (m=='lhs') then
    !multiply C^{1/2}
    alm(1,:,:) = clh*x
  else
    alm(1,:,:) = x
  end if

  !alm to map
  call alm2map(nside,lmax,lmax,alm,map)

  !multiply N^-1
  map = map*nij

  !map to alm
  call map2alm(nside,lmax,lmax,map,alm)

  !multiply C^{1/2}
  alm(1,:,:) = clh*alm(1,:,:)

  if (m=='lhs') then
    !make 1+C^{1/2}N^{-1}C^{1/2}
    v = x + alm(1,:,:)
  else
    v = alm(1,:,:)
  end if


end subroutine mat_multi

end module inpcurv


