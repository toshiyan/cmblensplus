!////////////////////////////////////////////////////!
! * Delensed BB
!////////////////////////////////////////////////////!

module delens
  use alkernel
  use delensing
  implicit none

contains


subroutine lintemplate(lmax,elmin,elmax,klmin,klmax,CE,Cm,WE,Wm,CB,gtype)
!*  Estimate of lensing template B-mode power spectrum (Wiener filters as inputs)
!*
!*  Args:
!*    :lmax (int)     : Maximum multipole of output spectrum
!*    :elmin (int)    : Minimum multipole of E
!*    :elmax (int)    : Maximum multipole of E
!*    :klmin (int)    : Minimum multipole of lensing mass
!*    :klmax (int)    : Maximum multipole of lensing mass
!*    :CE[l] (double) : Power spectrum of E-mode, with bounds (0:dlmax)
!*    :Cp[l] (double) : Power spectrum of lensing pontential, with bounds (0:dlmax)
!*    :WE[l] (double) : Wiener filter of E-mode, with bounds (0:dlmax)
!*    :Wp[l] (double) : Wiener filter of lensing potential, with bountd (0:dlmax)
!*
!*  Args(optional):
!*    :gtype (str) : specify type of the input Cp, p (default) or k.
!*
!*  Returns:
!*    :CB[l] (double) : Lensing B-mode power spectrum, with bounds (0:lmax)
!*
  implicit none
  !f2py intent(in) gtype, lmax, elmin, elmax, klmin, klmax, CE, WE, Cm, Wm
  !f2py intent(out) CB
  !f2py depend(elmax) CE, WE
  !f2py depend(klmax) Cm, Wm
  !f2py depend(lmax) CB
  !I/O
  character(*), intent(in) :: gtype
  integer, intent(in) :: lmax, elmin, elmax, klmin, klmax
  double precision, intent(in), dimension(0:elmax) :: CE, WE
  double precision, intent(in), dimension(0:klmax) :: Cm, Wm
  double precision, intent(out), dimension(0:lmax) :: CB
  !internal
  integer :: l
  double precision, dimension(0:klmax) :: Cp
  double precision, dimension(elmin:elmax) :: W1
  double precision, dimension(klmin:klmax) :: W2
  !opt4py :: gtype = 'p'

  Cp = 0d0
  if (gtype=='k') then
    do l = 1, klmax
      Cp(l) = 4d0*Cm(l)/dble(l**2+l)**2
    end do
  else
    Cp = Cm
  end if

  W1 = CE(elmin:elmax)*WE(elmin:elmax)
  W2 = Cp(klmin:klmax)*Wm(klmin:klmax)
  call conv_egrad((/1,lmax/),(/elmin,elmax/),(/klmin,klmax/),W1,W2,CB(1:lmax))

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
  !f2py intent(in) lmax, dlmin, dlmax, CE, Cp
  !f2py intent(out) CB
  !f2py depend(dlmax) CE, Cp
  !f2py depend(lmax) CB
  !I/O
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp
  double precision, intent(out), dimension(0:lmax) :: CB

  CB(0) = 0d0
  call conv_egrad((/1,lmax/),(/dlmin,dlmax/),(/dlmin,dlmax/),CE(dlmin:dlmax),Cp(dlmin:dlmax),CB(1:lmax))

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
  !f2py intent(in) lmax, dlmin, dlmax, CE, CB, Cp, NP1, NP2, Ag
  !f2py intent(out) DB
  !f2py depend(dlmax) CE, CB, Cp, NP1, NP2, Ag
  !f2py depend(lmax) DB
  !I/O
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE,CB,Cp,NP1,NP2,Ag
  double precision, intent(out), dimension(0:lmax) :: DB

  DB(0) = 0d0
  call Delensing_Bias((/1,lmax/),(/dlmin,dlmax/),CE(1:),CB(1:),Cp(1:),DB(1:),NP1(1:),Ag(1:),NP2(1:))

end subroutine delensbias_dom



end module delens


