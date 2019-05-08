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
!*    :lmax (int)     : Maximum multipole of residual ClBB
!*    :dlmin (int)    : Minimum multipole of E and lensing for delensing
!*    :dlmax (int)    : Maximum multipole of E and lensing for delensing
!*    :CE[l] (double) : Power spectrum of E-mode, with bounds (0:dlmax)
!*    :Cp[l] (double) : Power spectrum of lensing pontential, with bounds (0:dlmax)
!*    :WE[l] (double) : Wiener filter of E-mode, with bounds (0:dlmax)
!*    :Wp[l] (double) : Wiener filter of lensing potential, with bountd (0:dlmax)
!*
!*  Returns:
!*    :CB[l] (double) : Residual B-mode spectrum, with bounds (0:lmax)
!*
  implicit none
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp, WE, Wp
  double precision, intent(out), dimension(0:lmax) :: CB

  CB(0) = 0d0
  call res_clbb((/1,lmax/),(/dlmin,dlmax/),CE(1:dlmax),Cp(1:dlmax),CB(1:lmax),WE(1:dlmax),Wp(1:dlmax))

end subroutine resbb


subroutine lintemplate(lmax,dlmin,dlmax,CE,Cp,WE,Wp,CB)
!*  Estimate of lensing template B-mode power spectrum (Wiener filters as inputs)
!*
!*  Args:
!*    :lmax (int)     : Maximum multipole of output spectrum
!*    :dlmin (int)    : Minimum multipole of E and lensing for delensing
!*    :dlmax (int)    : Maximum multipole of E and lensing for delensing
!*    :CE[l] (double) : Power spectrum of E-mode, with bounds (0:dlmax)
!*    :Cp[l] (double) : Power spectrum of lensing pontential, with bounds (0:dlmax)
!*    :WE[l] (double) : Wiener filter of E-mode, with bounds (0:dlmax)
!*    :Wp[l] (double) : Wiener filter of lensing potential, with bountd (0:dlmax)
!*
!*  Returns:
!*    :CB[l] (double) : Lensing B-mode power spectrum, with bounds (0:lmax)
!*
  implicit none
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp, WE, Wp
  double precision, intent(out), dimension(0:lmax) :: CB
! [internal]
  double precision, dimension(dlmin:dlmax) :: W1, W2

  W1 = CE(dlmin:dlmax)*WE(dlmin:dlmax)
  W2 = Cp(dlmin:dlmax)*Wp(dlmin:dlmax)
  call conv_egrad((/1,lmax/),(/dlmin,dlmax/),W1,W2,CB)

end subroutine lintemplate


subroutine lensingbb(lmax,dlmin,dlmax,CE,Cp,CB)
!* Lensing B-mode power spectrum as a convolution of ClEE and Clpp
!*
!*  Args:
!*    :lmax (int)     : Maximum multipole of output spectrum
!*    :dlmin (int)    : Minimum multipole of E and lensing for delensing
!*    :dlmax (int)    : Maximum multipole of E and lensing for delensing
!*    :CE[l] (double) : Power spectrum of E-mode, with bounds (0:dlmax)
!*    :Cp[l] (double) : Power spectrum of lensing pontential, with bounds (0:dlmax)
!*
!*  Returns:
!*    :CB[l] (double) : Lensing B-mode power spectrum, with bounds (0:lmax)
!*
  implicit none
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp
  double precision, intent(out), dimension(0:lmax) :: CB

  CB(0) = 0d0
  call conv_egrad((/1,lmax/),(/dlmin,dlmax/),CE(dlmin:dlmax),Cp(dlmin:dlmax),CB(1:lmax))

end subroutine lensingbb


subroutine delensbias_dom(lmax,dlmin,dlmax,CE,CB,Cp,NP1,NP2,Ag,DB)
!*  Dominant term of the delensing bias in the B-mode internal delensing
!*
!*  Args:
!*    :lmax (int)     : Maximum multipole of output spectrum
!*    :dlmin (int)    : Minimum multipole of E and lensing for delensing
!*    :dlmax (int)    : Maximum multipole of E and lensing for delensing
!*    :CE[l] (double) : Power spectrum of E-mode, with bounds (0:dlmax)
!*    :CB[l] (double) : Power spectrum of B-mode, with bounds (0:dlmax)
!*    :Cp[l] (double) : Power spectrum of lensing pontential, with bounds (0:dlmax)
!*    :NP1[l] (double): Pol. noise spectrum for lensing reconstruction, with bounds (0:dlmax)
!*    :NP2[l] (double): Pol. noise spectrum for B-mode to be delensed, with bounds (0:dlmax)
!*    :Ag[l] (double) : Lensing reconstruction noise, with bounds (0:dlmax)
!*
!*  Returns:
!*    :DB[l] (double) : Lensing B-mode power spectrum, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE,CB,Cp,NP1,NP2,Ag
  double precision, intent(out), dimension(0:lmax) :: DB

  DB(0) = 0d0
  call Delensing_Bias((/1,lmax/),(/dlmin,dlmax/),CE(1:),CB(1:),Cp(1:),DB(1:),NP1(1:),Ag(1:),NP2(1:))

end subroutine delensbias_dom


end module delens

