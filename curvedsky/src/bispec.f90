!////////////////////////////////////////////////////!
! Bispectrum calcuation in fullsky
!////////////////////////////////////////////////////!

module bispec
  use alm_tools, only: alm2map, map2alm
  use constants, only: pi

  private alm2map
  private pi

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
  double complex, intent(in), dimension(1,0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(1,0:lmax,0:lmax) :: qlm
  !internal
  integer :: nside
  double precision, allocatable :: map(:)

  nside = int(lmax/2)

  allocate(map(0:12*nside**2-1))

  call alm2map(nside,lmax,lmax,alm,map)
  map = map + map**2
  call map2alm(nside,lmax,lmax,map,qlm)

  deallocate(map)

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
  character(4), intent(in), optional :: bstype
  integer, intent(in), optional :: bst
  integer, intent(in), optional, dimension(2) :: sL
  !f2py character(4) :: bstype = 'equi'
  !f2py integer :: bst = 2
  !f2py integer :: sL = 0
  !docstr :: sL = [int(bp[0]),int(bp[1])]
  !internal
  character(4) :: btype
  integer :: l1, l, lmax, b, bst0, aL(2), sL0(2), eL(2)
  double complex, allocatable :: klm(:,:)

  btype = bstype
  if (present(bstype))  btype = bstype

  bst0 = 2
  if (present(bst)) bst0 = bst

  !sque
  sL0  = int(bp(1:2))
  if (present(sL).and.sum(sL)/=0) sL0 = sL
  if (btype=='sque'.and..not.present(sL)) write(*,*) 'warning (bispec): sL is set to the lowest multipole bin'

  !isos
  aL = int(bp(bn/2:bn/2+1)) !middle bin

  !alm
  lmax = bp(bn+1)
  allocate(klm(0:lmax,0:lmax)); klm=0d0

  !if (present(alm).or.sum(abs(alm))/=0) then
  !  if (size(alm,dim=1)/=lmax+1.or.size(alm,dim=2)/=lmax+1) stop 'error (bispec): size of alm is not strange'
  !  klm = alm
  !else
    do l = 1, lmax
      klm(l,0) = dsqrt(2d0*l+1d0)
    end do
  !end if

  !compute binned bispectrum
  do b = 1, bn

    eL = int(bp(b:b+1))

    select case(btype)
    case('equi')
      call equi(eL(1),eL(2),klm(0:eL(2),0:eL(2)),norm(b),bst0)
    case('fold')
      call fold(eL(1),eL(2),klm(0:eL(2),0:eL(2)),norm(b),bst0)
    case('sque')
      l1 = max(eL(2),sL(2))
      call sque(eL,sL,l1,klm(0:l1,0:l1),norm(b),bst0)
    case('isos')
      l1 = max(eL(2),aL(2))
      call isos(eL,aL,l1,klm(0:l1,0:l1),norm(b),bst0)
    end select

  end do

  deallocate(klm)

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
  character(4), intent(in), optional :: bstype
  integer, intent(in), optional :: bst
  integer, intent(in), optional, dimension(2) :: sL
  !f2py integer :: bst = 2
  !f2py character :: bstype = 'equi'
  !f2py integer :: sL = 0
  !docstr :: sL = [int(bp[0]),int(bp[1])]
  !internal
  integer :: aL(2), b, eL(2), bst0, sL0(2), l1

  bst0 = 2
  if (present(bst)) bst0 = bst
  if (bp(bn+1)>lmax) stop 'error (equi_bin): not enough size of alm'

  sL0  = int(bp(1:2))
  if (present(sL).and.sum(sL)/=0) sL0 = sL

  aL = int(bp(bn/2:bn/2+1)) !middle bin

  do b = 1, bn
    eL = int(bp(b:b+1))
    select case (bstype)
    case ('equi')
      call equi(eL(1),eL(2),alm(0:eL(2),0:eL(2)),bis(b),bst0)
    case ('fold')
      call fold(eL(1),eL(2),alm(0:eL(2),0:eL(2)),bis(b),bst0)
    case ('sque')
      l1 = max(eL(2),sL(2))
      call sque(eL,sL0,l1,alm(0:l1,0:l1),bis(b),bst0)
    case ('isos')
      l1 = max(eL(2),aL(2))
      call isos(eL,aL,l1,alm(0:l1,0:l1),bis(b),bst0)
    end select
  end do

end subroutine bispec_bin


subroutine equi(lmin,lmax,alm,bispec,bst)
!*  Compute equilateral shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l,l]
!*
!*  Args:
!*    :lmin (int)        : Minimum multipole of the bin
!*    :lmax (int)        : Maximum multipole of the bin
!*    :alm [l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:lmax,0:lmax).
!*
!*  Args(optional):
!*    :bst (int)         : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*
!*  Returns:
!*    :bispec (double)   : Unnormalized binned reduced bispectrum at the bin, [lmin,lmax]
!*
  implicit none
  !I/O
  integer, intent(in) :: lmin, lmax
  integer, intent(in), optional :: bst
!f2py integer :: bst = 2
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:)
  double complex, allocatable :: klm(:,:,:)

  nside = 2*lmax
  if (present(bst))  nside = bst*lmax

  allocate(kmap(0:12*nside**2-1),klm(1,0:lmax,0:lmax))

  klm = 0d0
  do l = lmin, lmax !multipole filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,lmax,lmax,klm,kmap)
  bispec = sum(kmap**3)*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine equi


subroutine fold(lmin,lmax,alm,bispec,bst)
!*  Compute folded shape of the unnormalized binned reduced bispectrum for a given alm, b[l,l/2,l/2]
!*
!*  Args:
!*    :lmin (int)        : Minimum multipole of the bin
!*    :lmax (int)        : Maximum multipole of the bin
!*    :alm [l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:lmax,0:lmax).
!*
!*  Args(optional):
!*    :bst (int)         : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*
!*  Returns:
!*    :bispec (double)   : Unnormalized binned reduced bispectrum at the bin, [lmin,lmax]
!*
  implicit none
  !I/O
  integer, intent(in) :: lmin, lmax
  integer, intent(in), optional :: bst
!f2py integer :: bst = 2
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nside = 2*lmax
  if (present(bst))  nside = bst*lmax

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:lmax,0:lmax))

  klm = 0d0
  do l = lmin, lmax !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = max(2,int(lmin/2d0)), int(lmax/2d0)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,lmax,lmax,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,lmax,lmax,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine fold


subroutine sque(eL,sL,l1,alm,bispec,bst)
!*  Compute squeezed shape of the unnormalized binned reduced bispectrum for a given alm, b[sL,eL,eL]
!*
!*  Args:
!*    :eL[2] (int)        : Minimum and maximum multipoles of the bin, with bounds (2)
!*    :sL[2] (int)        : Minimum and maximum multipoles of the fixed bin, with bounds (2)
!*    :l1 (int)           : Maximum multipole of the input alm, satisfying eLmax,sLmax<=l1
!*    :alm [l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:lmax,0:lmax).
!*
!*  Args(optional):
!*    :bst (int)         : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*
!*  Returns:
!*    :bispec (double)   : Unnormalized binned reduced bispectrum at the bin, [lmin,lmax]
!*
  implicit none
  !I/O
  integer, intent(in) :: l1
  integer, intent(in), dimension(2) :: eL, sL
  integer, intent(in), optional :: bst
!f2py integer :: bst = 2
  double complex, intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  if (max(sL(2),eL(2))>l1) stop 'error (sque): l1 is too small'

  nside = 2*l1
  if (present(bst))  nside = bst*l1

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:l1,0:l1))

  klm = 0d0
  do l = sL(1), sL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = eL(1), eL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,l1,l1,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,l1,l1,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine sque


subroutine isos(eL,aL,l1,alm,bispec,bst)
!*  Compute isosceles shape of the unnormalized binned reduced bispectrum for a given alm, b[eL,aL,aL]
!*
!*  Args:
!*    :eL[2] (int)        : Minimum and maximum multipoles of the bin, with bounds (2)
!*    :aL[2] (int)        : Minimum and maximum multipoles of the fixed bin, with bounds (2)
!*    :l1 (int)           : Maximum multipole of the input alm, satisfying eLmax,sLmax<=l1
!*    :alm [l,m] (dcmplx) : Input harmonic coefficients, with bounds (0:lmax,0:lmax).
!*
!*  Args(optional):
!*    :bst (int)         : A parameter, bst=nside/lmax, which controls the accuracy of the calculation, default to 2. More accurate for a larger value.
!*
!*  Returns:
!*    :bispec (double)   : Unnormalized binned reduced bispectrum at the bin, [lmin,lmax]
!*
  implicit none
  !I/O
  integer, intent(in) :: l1
  integer, intent(in), dimension(2) :: eL, aL
  integer, intent(in), optional :: bst
!f2py integer :: bst = 2
  double complex, intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  if (max(aL(2),eL(2))>l1) stop 'error (isos): l1 is too small'

  nside = 2*l1
  if (present(bst))  nside = bst*l1

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:l1,0:l1))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = aL(1), aL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,l1,l1,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,l1,l1,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine isos


end module bispec


