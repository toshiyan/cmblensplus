!////////////////////////////////////////////////////!
! Bispectrum calcuation in fullsky
!////////////////////////////////////////////////////!

module bispec
  use alm_tools, only: alm2map, map2alm
  use constants, only: pi
  use hp_bsp,    only: equi, fold, sque, isos, xequi, xfold, xsque, xisos

  private alm2map
  private pi
  private equi, fold, sque, isos, xequi, xfold, xsque, xisos

contains 


subroutine make_quad_gauss(lmax,alm,qlm)
!*  Return a non-Gaussian alm. The corresponding non-Gaussian field is defined as
!*
!*    delta^NL(n) = delta^L(n) + delta^L(n)**2
!*
!*  where delta^L(n) is a gaussian field obtained from the input alm.
!*
!*  Args:
!*    :lmax (int)         : Maximum multipole of alm
!*    :alm [l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:lmax,0:lmax).
!*
!*  Returns:
!*    :qlm [l,m] (dcmplx) : Output harmonic coefficients of the non-Gaussian fields, with bounds (0:lmax,0:lmax).
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: qlm
  !internal
  integer :: nside
  double precision, allocatable :: map(:)
  double complex, allocatable :: alm_tmp(:,:,:)

  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))

  allocate(map(0:12*nside**2-1),alm_tmp(1,0:lmax,0:lmax))
  alm_tmp(1,:,:) = alm
  call alm2map(nside,lmax,lmax,alm_tmp,map)
  map = map + map**2
  call map2alm(nside,lmax,lmax,map,alm_tmp)
  qlm = alm_tmp(1,:,:)
  deallocate(map,alm_tmp)

end subroutine make_quad_gauss


subroutine bispec_norm(bn,bp,norm,bstype,bst,sL)
!*  Return normalization of the binned reduced bispectrum for a given multipole bin
!*
!*  Args:
!*    :bn (int)           : Number of multipole bins
!*    :bp [edge] (double) : Bin edges, with bounds (bn+1)
!*
!*  Args(optional):
!*    :bstype (str)  : Configuration of the bispectrum, default to equilateral
!*    :bst (int)     : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*    :sL[2] (int)   : The fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin
!*
!*  Returns:
!*    :norm [bin] (double) : Normalization of the binned reduced bispectrum at each bin, with bounds (bn)
!*
  implicit none
  !I/O
  integer, intent(in) :: bn
  double precision, intent(in), dimension(bn+1) :: bp
  double precision, intent(out), dimension(bn) :: norm
  !optional
  character(4), intent(in) :: bstype
  integer, intent(in) :: bst
  integer, intent(in), dimension(2) :: sL
  !opt4py :: bstype = 'equi'
  !opt4py :: bst = 2
  !opt4py :: sL = [0,0]
  !internal
  integer :: l, lmax, b, aL(2), sL0(2), eL(2)
  double precision :: nbisp
  double complex, allocatable :: klm(:,:)

  !sque
  sL0  = sL
  if (sum(sL)==0)  sL0 = int(bp(1:2))

  !isos
  aL = int(bp(bn/2:bn/2+1)) !middle bin

  !compute binned bispectrum
  do b = 1, bn

    eL = int(bp(b:b+1))

    !define maximum l
    select case(bstype)
    case ('equi')
      lmax = eL(2)
    case ('fold')
      lmax = eL(2)
    case ('sque')
      lmax = max(eL(2),sL(2))
    case ('isos')
      lmax = max(eL(2),aL(2))
    end select
  
    allocate(klm(0:lmax,0:lmax)); klm=0d0
    do l = 1, lmax
      klm(l,0) = dsqrt(2d0*l+1d0)
    end do

    select case(bstype)
    case('equi')
      call equi(eL(1),eL(2),klm,norm(b),bst)
    case('fold')
      call fold(eL(1),eL(2),klm,norm(b),bst)
    case('sque')
      call sque(eL,sL,lmax,klm,norm(b),bst)
    case('isos')
      call isos(eL,aL,lmax,klm,norm(b),bst)
    end select

    deallocate(klm)

  end do


end subroutine bispec_norm


subroutine bispec_bin(bn,bp,lmax,alm,bis,bstype,bst,sL)
!*  Return the unnormalized binned reduced bispectrum for a given multipole bin
!*
!*  Args:
!*    :bn (int)           : Number of multipole bins
!*    :bp [edge] (double) : Bin edges, with bounds (bn+1)
!*    :lmax (int)         : Maximum multipole of the input alm
!*    :alm [l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :bstype (str)  : Configuration of the bispectrum, default to equilateral
!*    :bst (int)     : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*    :sL[2] (int)   : The fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin
!*
!*  Returns:
!*    :bis [bin] (double) : The unnormalized binned reduced bispectrum at each bin, with bounds (bn)
!*
  implicit none
  !I/O
  integer, intent(in) :: bn, lmax
  double precision, intent(in), dimension(bn+1) :: bp
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double precision, intent(out), dimension(bn) :: bis
  !optional
  character(4), intent(in) :: bstype
  integer, intent(in) :: bst
  integer, intent(in), dimension(2) :: sL
  !opt4py :: bst = 2
  !opt4py :: bstype = 'equi'
  !opt4py :: sL = [0,0]
  !internal
  integer :: aL(2), b, eL(2), sL0(2), ilmax
  double complex, allocatable :: alm_tmp(:,:)

  if (bp(bn+1)>lmax) stop 'error (equi_bin): not enough size of alm'

  sL0  = sL
  if (sum(sL)==0)  sL0 = int(bp(1:2))

  aL = int(bp(bn/2:bn/2+1)) !middle bin

  do b = 1, bn
    eL = int(bp(b:b+1))
    
    !define maximum l
    select case(bstype)
    case ('equi','fold')
      ilmax = eL(2)
    case ('sque')
      ilmax = max(eL(2),sL(2))
    case ('isos')
      ilmax = max(eL(2),aL(2))
    end select
  
    allocate(alm_tmp(0:ilmax,0:ilmax))
    alm_tmp = alm(0:ilmax,0:ilmax)

    select case (bstype)
    case ('equi')
      call equi(eL(1),eL(2),alm_tmp,bis(b),bst)
    case ('fold')
      call fold(eL(1),eL(2),alm_tmp,bis(b),bst)
    case ('sque')
      call sque(eL,sL0,ilmax,alm,bis(b),bst)
    case ('isos')
      call isos(eL,aL,ilmax,alm,bis(b),bst)
    end select

    deallocate(alm_tmp)
  
  end do


end subroutine bispec_bin


subroutine xbispec_bin(bn,bp,lmax,n,alm,bis,bstype,bst,sL)
!*  Return the unnormalized binned reduced cross-bispectrum for a given multipole bin
!*
!*  Args:
!*    :bn (int)           : Number of multipole bins
!*    :bp [edge] (double) : Bin edges, with bounds (bn+1)
!*    :lmax (int)         : Maximum multipole of the input alm
!*    :n (int)            : Number of alms
!*    :alm [n,l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:n-1,0:lmax,0:lmax)
!*
!*  Args(optional):
!*    :bstype (str)  : Configuration of the bispectrum, default to equilateral
!*    :bst (int)     : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*    :sL[2] (int)   : The fixed bin for the squeezed configuration, b[sL,eL,eL], default to the lowest multipole bin
!*
!*  Returns:
!*    :bis [bin] (double) : The unnormalized binned reduced bispectrum at each bin, with bounds (bn)
!*
  implicit none
  !I/O
  character(4), intent(in) :: bstype
  integer, intent(in) :: bn, lmax, bst, n
  integer, intent(in), dimension(2) :: sL
  double precision, intent(in), dimension(bn+1) :: bp
  double complex, intent(in), dimension(1:n,0:lmax,0:lmax) :: alm
  double precision, intent(out), dimension(bn) :: bis
  !opt4py :: bst = 2
  !opt4py :: bstype = 'equi'
  !opt4py :: sL = [0,0]
  !internal
  integer :: aL(2), b, eL(2), l1, sL0(2)

  if (bp(bn+1)>lmax) stop 'error (xbispec_bin): not enough size of alm'

  sL0  = sL
  if (sum(sL)==0)  sL0 = int(bp(1:2))

  aL = int(bp(bn/2:bn/2+1)) !middle bin

  do b = 1, bn
    eL = int(bp(b:b+1))
    select case (bstype)
    case ('equi')
      call xequi(eL(1),eL(2),n,alm(:,0:eL(2),0:eL(2)),bis(b),bst)
    case ('fold')
      call xfold(eL(1),eL(2),n,alm(:,0:eL(2),0:eL(2)),bis(b),bst)
    case ('sque')
      l1 = max(eL(2),sL0(2))
      call xsque(eL,sL0,l1,n,alm(:,0:l1,0:l1),bis(b),bst)
    case ('isos')
      l1 = max(eL(2),aL(2))
      call xisos(eL,aL,l1,n,alm(:,0:l1,0:l1),bis(b),bst)
    end select
  end do

end subroutine xbispec_bin


end module bispec


