!////////////////////////////////////////////////////!
! * LAPACK
!////////////////////////////////////////////////////!


module lapack
  implicit none


contains


subroutine solve_eigen_problem(n,W,CE,CB)
  implicit none
  !I/O
  integer, intent(in) :: n
  double precision, intent(out) :: W(1:n)
  double precision, intent(inout) :: CE(1:n,1:n), CB(1:n,1:n)
  !internal
  integer :: lwork, info, i, c1, c2
  !double precision :: rate, t1, t2
  double precision, allocatable :: work(:)

  lwork = 3*n-1
  allocate(work(1:lwork))
  !call system_clock(c1,rate)
  !call cpu_time(t1)
  CALL DSYGV(1,'V','U',n,CB,n,CE,n,W,work,lwork,info)
  !call system_clock(c2,rate)
  !call cpu_time(t2)
  !write(*,*) "real time (eigen solver)=", (c2-c1)/rate
  !write(*,*) "cpu time (eigen solver)=", t2-t1
  deallocate(work)

  open(unit=20,file="lambda.dat",status="replace")
  do i = 1, n
    write(20,'(I6,1X,(E12.5))') i, W(i)
  end do
  close(20)

  if(.not.info==0) then
    write(*,*) "info =", info
    stop "stop with error in solving generalized eigenvalue problem"
  end if

end subroutine solve_eigen_problem


end module lapack

