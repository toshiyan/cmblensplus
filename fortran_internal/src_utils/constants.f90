!///////////////////////////////////////////////////////!
! * Constants
!///////////////////////////////////////////////////////!

module constants
  implicit none

  !* precision
  integer, parameter :: dl  = KIND(1d0)
  integer, parameter :: dlc = KIND(0d0)

  !* constants
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0
  double precision, parameter :: ac2rad = pi/180d0/60d0
  double precision, parameter :: twopi=2d0*pi, fourpi=4d0*pi
  double precision, parameter :: const_c = 2.99792458d8
  double precision, parameter :: const_c_km = 2.99792458d5
  double precision, parameter :: const_hP = 6.62606896d-34
  double precision, parameter :: const_G = 6.67428d-11
  double precision, parameter :: const_Gyr = 3.1556926d16
  double precision, parameter :: Tcmb = 2.726d6  ! micro K
  complex(dlc), parameter :: iu = (0d0,1d0)
  complex(dlc), parameter :: i1 = (1d0,0d0)

  !* statistical significance
  double precision, parameter :: signif(3) = [0.68269d0,0.95450d0,0.99730d0]

  !* CMB spectra
  integer, parameter :: TT=1, TE=2, EE=3, BB=4, dd=5, Td=6, Ed=7, oo=8

end module constants

