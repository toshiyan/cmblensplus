!////////////////////////////////////////////////////!
! * Delensed BB
!////////////////////////////////////////////////////!

module delens
  use alkernel
  use rec_lens, only: qeb
  use delensing
  implicit none

contains


subroutine qeb_iter(rlmin,rlmax,dlmin,dlmax,fCEE,Cpp,oEE,oBB,Alg,Alc,iter,conv)
  implicit none
  !I/O
  integer, intent(in) :: rlmin, rlmax, dlmin, dlmax
  double precision, intent(in), dimension(0:rlmax) :: fCEE, oEE, oBB
  double precision, intent(in), dimension(0:dlmax) :: Cpp
  double precision, intent(out), dimension(0:dlmax) :: Alg, Alc
  !optional
  integer, intent(in), optional :: iter
  double precision, intent(in), optional :: conv
  !f2py integer :: iter = 1
  !f2py double precision :: conv = 0.001
  !internal
  integer :: i, n, l
  double precision :: ratio
  double precision, dimension(:), allocatable :: AlgEB, rCBB

  allocate(AlgEB(0:dlmax),rCBB(0:rlmax));  AlgEB=0d0;  rCBB=0d0

  if (rlmax<dlmax) stop 'error (qeb_iter): does not support rlmax<dlmax case'

  !initial values
  ratio = 1d0
  rCBB  = oBB(:rlmax)

  do n = 1, iter !loop for iteration 

    !* lensing reconstruction with EB
    call qeb(dlmax,rlmin,rlmax,fCEE,oEE,rCBB,AlgEB,Alc)

    !* convergence check
    if (n>=2) then
      ratio = (sum(Alg)/sum(AlgEB)-1d0)/dble(dlmax)
      write(*,*) n, ratio
    end if
    Alg = AlgEB

    if (abs(ratio) < conv) exit

    !* delensing with EB-estimator
    call CLBB_EST((/rlmin,rlmax/),(/dlmin,dlmax/),fCEE(1:dlmax),Cpp(1:dlmax),oEE(1:dlmax)-fCEE(1:dlmax),AlgEB,rCBB)
    rCBB = oBB(:rlmax) - rCBB !delensed B-mode

    if(n==iter) stop 'not converged'

  end do

  deallocate(AlgEB,rCBB)

end subroutine qeb_iter


subroutine resbb(lmax,dlmin,dlmax,CE,Cp,WE,Wp,CB)
!* residual ClBB = ClBB^lin - ClBB^est
!
  implicit none
! [inputs]  
!   lmax --- maximum multipole of residual ClBB
!   dL --- multipole range of delensing
!   CE, Cp --- power spectrum of E-mode and lensing pontential
!   WE, Wp --- Wiener filters of E-mode and lensing potential
  integer, intent(in) :: lmax, dlmin, dlmax
  double precision, intent(in), dimension(0:dlmax) :: CE, Cp, WE, Wp
!
! [outputs]
!   CB --- residual B-mode
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

