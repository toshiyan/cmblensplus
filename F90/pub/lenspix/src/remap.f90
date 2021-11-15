
module remap_cmb
  use spinalm_tools
  implicit none

contains


subroutine simple_remapping(nside,lmax,alm,grad,tqu)
  implicit none
  integer, intent(in) :: nside, lmax
  complex, intent(in) :: alm(:,:,:), grad(:)
  real, intent(out)   :: tqu(:,:)
  !internal
  Type(HealpixInfo)  :: H
  real :: interp_factor

  interp_factor = 1.5*(2048.0/real(nside))
  if (nside<1024)  interp_factor=1.5

  call HealpixInit(H,nside,lmax,.true.,w8dir='',method=division_equalrows) 
  call alm2LensedmapInterpCyl(H,lmax,alm,grad,TQU,interp_factor)
  call HealpixFree(H)

end subroutine simple_remapping


end module remap_cmb

