!////////////////////////////////////////////////////!
! Interface of Spherical Harmonic Transform
!////////////////////////////////////////////////////!

module hp_spht
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

  tmap = 0d0
  select case(n)
  case(1)
    call alm2map(nside,lmax,mmax,alm,tmap(0:npix-1,1))
  case(2)
    call alm2map_spin(nside,lmax,mmax,2,alm(1:2,0:lmax,0:mmax),tmap(0:npix-1,1:2))
  case(3)
    call alm2map(nside,lmax,mmax,alm(1:1,0:lmax,0:mmax),tmap(0:npix-1,1))
    call alm2map_spin(nside,lmax,mmax,2,alm(2:3,0:lmax,0:mmax),tmap(0:npix-1,2:3))
  case default
    stop 'error (spht_alm2map): n should be 1, 2, or 3'
  end select

  do ni = 1, n
    map(ni,0:npix-1) = tmap(0:npix-1,ni)
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

  alm = 0d0
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


end module hp_spht


