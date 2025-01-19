
subroutine EST_BHE(est,AL,els,eL) 
! Transform biased estimator (est) to unbiased estimator (x)
!   est = R * x 
! where R is the response function. 
  implicit none 
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: AL(:,:), els(:)
  double complex, intent(inout) :: est(:,:)
  !internal
  integer :: n, m, i, j
  double precision, allocatable :: R(:,:)

  M = size(est,dim=1)
  allocate(R(M,M))

  do n = 1, size(est,dim=2)
    if(ran(els(n),eL)) then
      !* response function (not symmetric)
      do i = 1, M
        R(i,i) = 1d0
        do j = i + 1, M
          R(i,j) = Al(mdiag(i,M),n)*Al(mlab(i,j,M),n)
          R(j,i) = Al(mdiag(j,M),n)*Al(mlab(i,j,M),n)
        end do
      end do
      CALL INV_LAPACK(R)
      est(:,n) = matmul(R,est(:,n))
    end if
  end do

  deallocate(R)

end subroutine EST_BHE

