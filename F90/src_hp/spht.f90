!////////////////////////////////////////////////////!
! Spherical Harmonic Transform
!////////////////////////////////////////////////////!

module spht
  !from Healpix
  use alm_tools, only: alm2map, alm2map_spin, map2alm, map2alm_spin
  implicit none

  private alm2map, alm2map_spin, map2alm, map2alm_spin

contains


subroutine spht_alm2map(n,npix,lmax,mmax,alm,map)
  implicit none
  integer, intent(in) :: n, npix, lmax, mmax
  double complex, intent(in), dimension(n,0:lmax,0:mmax) :: alm
  double precision, intent(out), dimension(n,0:npix-1) :: map
  integer :: nside, ni
  double precision :: tmap(0:npix-1,1:n)

  nside = int(sqrt(npix/12d0))

  select case(n)
  case(1)
    call alm2map(nside,lmax,mmax,alm,tmap(:,1))
  case(2)
    call alm2map_spin(nside,lmax,mmax,2,alm(1:2,:,:),tmap(:,1:2))
  case(3)
    call alm2map(nside,lmax,mmax,alm(1:1,:,:),tmap(:,1))
    call alm2map_spin(nside,lmax,mmax,2,alm(2:3,:,:),tmap(:,2:3))
  case default
    stop 'error (spht_alm2map): n should be 1, 2, or 3'
  end select

  do ni = 1, n
    map(ni,:) = tmap(:,ni)
  end do

end subroutine spht_alm2map


subroutine spht_map2alm(n,npix,lmax,mmax,map,alm)
  implicit none
  integer, intent(in) :: n, npix, lmax, mmax
  double precision, intent(in), dimension(n,0:npix-1) :: map
  double complex, intent(out), dimension(n,0:lmax,0:mmax) :: alm
  integer :: nside, ni
  double precision :: tmap(0:npix-1,1:n)

  nside = int(sqrt(npix/12d0))

  do ni = 1, n
    tmap(:,ni) = map(ni,:)
  end do

  select case(n)
  case(1)
    call map2alm(nside,lmax,mmax,tmap(:,1),alm)
  case(2)
    call map2alm_spin(nside,lmax,mmax,2,tmap(:,1:2),alm(1:2,:,:))
  case(3)
    call map2alm(nside,lmax,mmax,tmap(:,1),alm(1:1,:,:))
    call map2alm_spin(nside,lmax,mmax,2,tmap(:,2:3),alm(2:3,:,:))
  case default
    stop 'error (spht_map2alm): n should be 1, 2, or 3'
  end select

end subroutine spht_map2alm


subroutine mat_multi(n,npix,lmax,clh,nij,x,v,mtype)
!* multiplying matrix, 1 + C^1/2 N^-1 C^1/2, if mtype = lhs and C^1/2 N^-1 otherwise.
  implicit none
  !I/O
  character(3), intent(in) :: mtype
  integer, intent(in) :: n, npix, lmax
  double precision, intent(in), dimension(n,0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: x
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  integer :: p
  double precision :: map(n,0:npix-1)
  double complex :: alm(n,0:lmax,0:lmax)

  if (mtype=='lhs') then
    !multiply C^{1/2}
    alm = clh*x
  else
    alm = x
  end if

  !multiply noise covariance in pixel space
  call spht_alm2map(n,npix,lmax,lmax,alm,map)
  map = map*nij
  call spht_map2alm(n,npix,lmax,lmax,map,alm)

  !multiply C^{1/2}
  alm = clh*alm

  if (mtype=='lhs') then
    !make 1+C^{1/2}N^{-1}C^{1/2}
    v = x + alm
  else
    v = alm
  end if

end subroutine mat_multi


subroutine mat_multi_freq_rhs(n,k,npix,lmax,clh,nij,x,v)
  implicit none
  !I/O
  integer, intent(in) :: n, k, npix, lmax
  double precision, intent(in), dimension(n,k,0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(n,k,0:npix-1) :: nij, x
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  integer :: ki
  double precision :: map(n,k,0:npix-1)
  double complex :: alm(n,0:lmax,0:lmax)

  v   = 0d0
  map = x*nij
  do ki = 1, k
    call spht_map2alm(n,npix,lmax,lmax,map(:,ki,:),alm)
    v = v + clh(:,ki,:,:)*alm
  end do

end subroutine mat_multi_freq_rhs


subroutine mat_multi_freq_lhs(n,k,npix,lmax,clh,nij,x,v)
  implicit none
  !I/O
  integer, intent(in) :: n, k, npix, lmax
  double precision, intent(in), dimension(n,k,0:lmax,0:lmax) :: clh
  double precision, intent(in), dimension(n,k,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: x
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: v
  !internal
  integer :: p, ki
  double precision :: map(n,0:npix-1)
  double complex :: alm(n,0:lmax,0:lmax), blm(n,0:lmax,0:lmax)

  blm = 0d0
  do ki = 1, k
    alm = clh(:,ki,:,:)*x
    !multiply noise covariance in pixel space
    call spht_alm2map(n,npix,lmax,lmax,alm,map)
    map = map*nij(:,ki,:)
    call spht_map2alm(n,npix,lmax,lmax,map,alm)
    !multiply C^{1/2}
    blm = blm + clh(:,ki,:,:)*alm
  end do

  !make 1+C^{1/2}N^{-1}C^{1/2}
  v = x + blm

end subroutine mat_multi_freq_lhs


end module spht


