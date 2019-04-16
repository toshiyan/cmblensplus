!////////////////////////////////////////////////////!
! * LAPACK and LAPACK95 
!////////////////////////////////////////////////////!

module lapack95
  use f95_LAPACK
  implicit none

contains

subroutine INV_LAPACK(K,M,accuracy)
!* interface of the matrix inversion
  implicit none
  !I/O
  double precision, intent(inout) :: K(:,:)
  double precision, intent(in), optional :: accuracy
  double precision, intent(out), optional :: M(:,:)
  !internal
  integer :: N
  integer, allocatable :: IPIV(:)
  double precision, allocatable :: iK(:,:)

  N = size(K,dim=1)
  allocate(iK(N,N),IPIV(N))
  !* inverting matrix
  if(N==1) then
    iK = 1.d0/K
  else
    iK = K
    call LA_GETRF(iK,IPIV)
    call LA_GETRI(iK,IPIV)
  end if
  if(present(accuracy)) call CHECK_INV_LAPACK(K,iK,accuracy)
  if(present(M)) then
    M = iK
  else
    K = iK
  end if
  deallocate(iK,IPIV)

end subroutine INV_LAPACK


subroutine CHECK_INV_LAPACK(M,iM,accuracy)
!* input matrix (M) and its inverse (iM)
  implicit none
  !I/O
  double precision, intent(in) :: M(:,:), iM(:,:), accuracy
  !internal
  integer :: i, j, n
  double precision, allocatable :: lam(:), A(:,:), tA(:,:)

  n = size(M,dim=1)
  allocate(A(n,n),tA(n,n),lam(n))
  A = matmul(iM,M)  ! should be identity matrix
  do i = 1, n
    ! diagonal
    if (abs(A(i,i)-1.d0) > accuracy) then
      write(*,*) A(i,i)
      stop "matrix inversion error (diag)"
    end if
    do j = 1, n
      if (abs(A(i,j)) > accuracy .and.(.not.i==j)) then
        write(*,*) A(i,j)
        stop "matrix inversion error (off-diag)"
      end if
    end do
  end do
  deallocate(A,tA,lam)

end subroutine CHECK_INV_LAPACK


subroutine SYGV_LAPACK(K,U,D,M)
!* solve generalized eigenvalue problem: K*x = D*M*x
!* output ith eigenvector U(:,i) with ith eigenvalue D(i)
  implicit none
  !I/O
  double precision, intent(inout) :: K(:,:)
  double precision, intent(out) :: D(:), U(:,:)
  double precision, optional :: M(:,:)
  !internal
  integer :: i, j, n
  double precision, allocatable :: B(:,:)
 
  n = size(D)
  allocate(B(n,n))
  if(present(M)) then
    B = M
  else
    do i = 1, N
      do j = 1, N
        if (i==j) then
          B(i,j) = 1.d0
        else
          B(i,j) = 0
        end if
      end do
    end do
  end if
  U = K
  call LA_SYGV(U,B,D,1,'V')
  deallocate(B)
 
end subroutine SYGV_LAPACK



subroutine CHECK_SYGV_LAPACK(D,U,M)
  implicit none
  double precision, intent(in) :: D(:), U(:,:), M(:,:)
  double precision, allocatable :: C(:,:)
  integer :: n, i, j

  n = size(M,dim=1)
  allocate(C(n,n))

  C = matmul(transpose(U),matmul(M,U))

  do i = 1, n
    if(abs(D(i)-C(i,i))>1e-7) then
      write(*,*) i, D(i)-C(i,i)
      stop 'error in LU decomposition : eigenmode 1'
    end if
  end do

  do i = 1, n
    do j = i+1, n
      if(abs(C(i,j))>1e-7.and..not.i==j) then
        write(*,*) i, j, C(i,j)
        stop 'error in LU decomposition : eigenmode 2'
      end if
    end do
  end do

  deallocate(C)

end subroutine CHECK_SYGV_LAPACK


end module lapack95

