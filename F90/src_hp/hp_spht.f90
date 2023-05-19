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
  double precision, allocatable :: tmap(:), pmap(:,:)
  double complex, allocatable :: tlm(:,:,:), plm(:,:,:)

  nside = int(sqrt(npix/12d0))
  
  allocate(tmap(0:npix-1),pmap(0:npix-1,2),tlm(1,0:lmax,0:lmax),plm(2,0:lmax,0:lmax))

  tmap = 0d0
  select case(n)
  case(1)
    call alm2map(nside,lmax,mmax,alm,tmap)
    map(1,:) = tmap
  case(2)
    call alm2map_spin(nside,lmax,mmax,2,alm,pmap)
    map(1,:) = pmap(:,1)
    map(2,:) = pmap(:,2)
  case(3)
    tlm(1,:,:) = alm(1,:,:)
    call alm2map(nside,lmax,mmax,tlm,tmap)
    map(1,:) = tmap
    plm(1,:,:) = alm(2,:,:)
    plm(2,:,:) = alm(3,:,:)
    call alm2map_spin(nside,lmax,mmax,2,plm,pmap)
    map(2,:) = pmap(:,1)
    map(3,:) = pmap(:,2)
  case default
    stop 'error (spht_alm2map): n should be 1, 2, or 3'
  end select


end subroutine spht_alm2map


subroutine spht_map2alm(n,npix,lmax,mmax,map,alm)
  implicit none
  integer, intent(in) :: n, npix, lmax, mmax
  double precision, intent(in), dimension(n,0:npix-1) :: map
  double complex, intent(out), dimension(n,0:lmax,0:mmax) :: alm
  integer :: nside, ni
  double precision, allocatable :: tmap(:), pmap(:,:)
  double complex, allocatable :: tlm(:,:,:), plm(:,:,:)

  nside = int(sqrt(npix/12d0))

  allocate(tmap(0:npix-1),pmap(0:npix-1,2),tlm(1,0:lmax,0:lmax),plm(2,0:lmax,0:lmax))
  
  select case(n)
  case(1)
    tmap = map(1,:)
    call map2alm(nside,lmax,mmax,tmap,tlm)
    alm = tlm
  case(2)
    pmap(:,1) = map(1,:)
    pmap(:,2) = map(2,:)
    call map2alm_spin(nside,lmax,mmax,2,pmap,plm)
    alm = plm
  case(3)
    tmap = map(1,:)
    pmap(:,1) = map(2,:)
    pmap(:,2) = map(3,:)
    call map2alm(nside,lmax,mmax,tmap,tlm)
    call map2alm_spin(nside,lmax,mmax,2,pmap,plm)
    alm(1:1,:,:) = tlm
    alm(2:3,:,:) = plm
  case default
    stop 'error (spht_map2alm): n should be 1, 2, or 3'
  end select
  
  deallocate(tmap,pmap,tlm,plm)

end subroutine spht_map2alm


end module hp_spht


