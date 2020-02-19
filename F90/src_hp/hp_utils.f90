module hp_utils
  use general, only: check_error, str
  implicit none

  private check_error, str

contains


subroutine trans_alm2array_1d(n,lmax,arrn,alm,arr)
  implicit none
  integer, intent(in) :: n, lmax, arrn
  double complex, intent(in) :: alm(1:n,0:lmax,0:lmax)
  double complex, intent(out) :: arr(1:arrn)
  integer :: i, l, m, ni

  arr = 0d0
  i = 0
  do ni = 1, n
    do l = 0, lmax
      do m = 0, l
        i = i + 1
        call check_error(i>arrn,'wrong size')
        arr(i) = alm(ni,l,m)
      end do
    end do
  end do

end subroutine trans_alm2array_1d


subroutine trans_array2alm_1d(n,arrn,lmax,arr,alm)
  implicit none
  integer, intent(in) :: n, lmax, arrn
  double complex, intent(in)  :: arr(1:arrn)
  double complex, intent(out) :: alm(1:n,0:lmax,0:lmax)
  integer :: i, l, m, ni

  alm = 0d0
  i = 0
  do ni = 1, n
    do l = 0, lmax
      do m = 0, l
        i = i + 1
        call check_error(i>arrn,'wrong size (F90/src_hp/hp_utils.f90.trans_array2alm_1d)','i='//str(i)//',arrn='//str(arrn))
        alm(ni,l,m) = arr(i)
      end do
    end do
  end do

end subroutine trans_array2alm_1d


subroutine trans_array2alm_2d(n,lmax,arrn,arr,alm)
  implicit none
  integer, intent(in) :: n, lmax, arrn
  double complex, intent(in)  :: arr(1:arrn,1:arrn)
  double complex, intent(out) :: alm(1:n,0:lmax,0:lmax,1:n,0:lmax,0:lmax)
  integer :: i, j, l, m, p, q, ni, nj

  alm = 0d0
  i = 0
  do ni = 1, n
    do l = 0, lmax
      do m = 0, l
        i = i + 1
        call check_error(i>arrn,'wrong size (F90/src_hp/hp_utils.f90.trans_array2alm_2d)')
        j = 0
        do nj = 1, n
          do p = 0, lmax
            do q = 0, p
              j = j + 1
              call check_error(j>arrn,'wrong size (F90/src_hp/hp_utils.f90.trans_array2alm_2d)')
              alm(ni,l,m,nj,p,q) = arr(i,j)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine trans_array2alm_2d


end module hp_utils

