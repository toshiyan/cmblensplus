!////////////////////////////////////////////////////!
! * Delensed BB
!////////////////////////////////////////////////////!

module delens
  use alkernel
  use delensing
  implicit none

contains


subroutine resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp,CB)
!*  Residual B-mode spectrum; ClBB = ClBB^lin - ClBB^est
!*
!*  Args:
!*    lmax (int)     : maximum multipole of residual ClBB
!*    dlmin (int)    : minimum multipole of E and lensing for delensing
!*    dlmax (int)    : maximum multipole of E and lensing for delensing
!*    CE[l] (double) : power spectrum of E-mode, with bounds (0:dlmax)
!*    Cp[l] (double) : power spectrum of lensing pontential, with bounds (0:dlmax)
!*    WE[l] (double) : Wiener filter of E-mode, with bounds (0:dlmax)
!*    Wp[l] (double) : Wiener filter of lensing potential, with bountd (0:dlmax)
!*
!*  Returns:
!*    CB[l] (double) : residual B-mode spectrum, with bounds (0:lmax)
!*
  implicit none
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp, WE, Wp
  double precision, intent(out) :: CB(0:lmax)

  CB(0) = 0d0
  call res_clbb((/1,lmax/),(/dlmin,dlmax/),CE(1:dlmax),Cp(1:dlmax),CB(1:lmax),WE(1:dlmax),Wp(1:dlmax))

end subroutine resbb


subroutine lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp,Cl)
! * Estimate of lensing template B-mode power spectrum (Wiener filters as inputs)
  implicit none
!
! [input]
!   lmax --- maximum multipole of residual ClBB
!   dL --- multipole range of delensing
!   CE --- E-mode power spectrum
!   Cp --- lensing potential power spectrum
!   WE --- E-mode Wiener filter
!   Wp --- lensing potential Wiener filter
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp, WE, Wp
!
! [output]
!   Cl --- estimated lensing B-mode power spectrum
  double precision, intent(out) :: Cl(0:lmax)
!
! [internal]
  double precision, dimension(dlmin:dlmax) :: W1, W2

  W1 = CE(dlmin:dlmax)*WE(dlmin:dlmax)
  W2 = Cp(dlmin:dlmax)*Wp(dlmin:dlmax)
  call conv_egrad((/1,lmax/),(/dlmin,dlmax/),W1,W2,Cl)

end subroutine lintemplate


subroutine lensingbb(lmax,dlmin,dlmax,CE,Cp,CB)
! * Lensing B-mode power spectrum as a convolution of ClEE and Clpp
  implicit none
!
! [input]
!   lmax --- maximum multipole of residual ClBB
!   dL --- multipole range of convolution
!   CE, Cp --- power spectrum of E-mode and lensing pontential
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp
!
! [output]
!   CB --- residual B-mode
  double precision, intent(out) :: CB(0:lmax)

  CB(0) = 0d0
  call conv_egrad((/1,lmax/),(/dlmin,dlmax/),CE(dlmin:dlmax),Cp(dlmin:dlmax),CB(1:lmax))

end subroutine lensingbb


subroutine delensbias_dom(lmax,dlmin,dlmax,EE,BB,pp,NP1,NP2,Ag,DBl)
  implicit none
  !I/O
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: EE,BB,pp,NP1,NP2,Ag
  double precision, intent(out) :: DBL(0:lmax)

  DBL(0) = 0d0
  call Delensing_Bias((/1,lmax/),(/dlmin,dlmax/),EE(1:),BB(1:),pp(1:),DBL(1:),NP1(1:),Ag(1:),NP2(1:))

end subroutine delensbias_dom


end module delens

