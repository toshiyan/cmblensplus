subroutine alm2cl_spc(lmax,cl,alm1,alm2,norm)
  implicit none
  !I/O
  integer, intent(in) :: lmax
  complex, intent(in), dimension(0:lmax,0:lmax) :: alm1
  double precision, intent(out), dimension(0:lmax) :: cl
  !(optional)
  complex, intent(in), dimension(0:lmax,0:lmax), optional :: alm2
  double precision, intent(in), optional :: norm
  !f2py complex :: alm2 = 0
  !f2py double precision :: norm = 1
  !intenral
  integer :: l

  Cl = 0d0
  if (present(alm2).and..not.sum(abs(alm2))/=0) then
    do l = 1, lmax
      Cl(l) = ( dble(alm1(l,0)*alm2(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm2(l,1:l))))/(2.*l+1.)
    end do
  else
    do l = 1, lmax
      Cl(l) = ( dble(alm1(l,0)*alm1(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm1(l,1:l))))/(2.*l+1.)
    end do
  end if
  if (present(norm))  Cl = Cl*norm

end subroutine alm2cl_spc



