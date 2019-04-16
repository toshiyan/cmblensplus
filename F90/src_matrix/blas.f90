!////////////////////////////////////////////////////!
! * BLAS Call
!////////////////////////////////////////////////////!

module blas
  implicit none
  integer, parameter :: dlc = KIND(0d0)

  !interface operator (.KP.) !Kronecker product
  interface KP !Kronecker product
    module procedure Kronecker_Dot_CMPLX, Kronecker_Dot_DBLE
  end interface

  private dlc

contains


function BLAS_DGEMM(M1,M2)  result(MM)
! return matrix multiplication
!   M_{m,n} = sum_k M1_{m,k}M_{k,n}
  implicit none
  !I/O
  double precision, intent(in) :: M1(:,:), M2(:,:)
  double precision :: MM(size(M1,dim=1),size(M2,dim=2))
  !internal
  integer :: m,k,n

  m = size(M1,dim=1)
  k = size(M1,dim=2)
  n = size(M2,dim=2)
  if (.not.size(M2)==k) stop "can not perform matrix multiplication"
  call dgemm('N','N',m,n,k,1.d0,M1,m,M2,k,0.d0,MM,m)

end function BLAS_DGEMM


function BLAS_CGEMM(M1,M2)  result(MM)
! return matrix multiplication
!   M_{m,n} = sum_k M1_{m,k}M_{k,n}
  implicit none
  !I/O
  complex(dlc), intent(in) :: M1(:,:), M2(:,:)
  complex(dlc) :: MM(size(M1,dim=1),size(M2,dim=2))
  !internal
  integer :: m,k,n

  m = size(M1,dim=1)
  k = size(M1,dim=2)
  n = size(M2,dim=2)
  if (.not.size(M2)==k) stop "can not perform matrix multiplication"
  call dgemm('N','N',m,n,k,1.d0,M1,m,M2,k,0.d0,MM,m)

end function BLAS_CGEMM


function Kronecker_Dot_DBLE(a,b) result(M)
! return Kronecker product: M = (a)x(b)
  implicit none
  !I/O
  double precision, intent(in) :: a(:), b(:)
  double precision :: M(size(a),size(b))
  !internal
  double precision :: va(1:size(a),1:1), vb(1:1,1:size(b))

  va(:,1) = a
  vb(1,:) = b
  M = BLAS_DGEMM(va,vb)

end function Kronecker_Dot_DBLE


function Kronecker_Dot_CMPLX(a,b) result(M)
  implicit none
  complex(dlc), intent(in) :: a(:), b(:)
  complex(dlc) :: M(size(a),size(b))
  !internal
  complex(dlc) :: va(1:size(a),1:1), vb(1:1,1:size(b))

  va(:,1) = a
  vb(1,:) = b
  M = BLAS_CGEMM(va,vb)

end function Kronecker_Dot_CMPLX


end module blas

