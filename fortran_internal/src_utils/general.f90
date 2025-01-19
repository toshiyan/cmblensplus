!////////////////////////////////////////////////////!
! * Intrinsic Subtourintes
!////////////////////////////////////////////////////!

module general
  implicit none
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0
  double complex, parameter :: iu = (0d0,1d0)

  !* Gauss-Legendre Quadrature
  type gauss_legendre_params
    integer :: n
    double precision :: eps
    double precision, allocatable :: z(:), w(:)
  end type gauss_legendre_params

  interface filecolumns
    module procedure filecolumns_fname, filecolumns_unit
  end interface

  interface filelines
    module procedure filelines_fname, filelines_unit
  end interface

  interface loadtxt
    module procedure loadtxt_1d_d, loadtxt_1d_c, loadtxt_2d_d, loadtxt_2d_c
  end interface loadtxt

  interface savetxt
    module procedure savetxt_1d_d, savetxt_1d_c, savetxt_2d_d, savetxt_2d_c
  end interface savetxt

  interface str
    module procedure str_n, str_fix, str_fix_r, str_fix_d
  end interface

  interface linspace
    module procedure linspace_dble, linspace_int, linspace_int_simple, linspace_dble_div, linspace_int_div
  end interface

  interface interp_lin
    module procedure interp_lin_arr, interp_lin_narr
  end interface

  interface ave
    module procedure ave_1d, ave_2d
  end interface ave

  interface save_average
    module procedure save_average_2d
  end interface save_average

  interface vec2mat
    module procedure vec2mat_dble, vec2mat_cmplx, vec2mat_diag
  end interface 

  private pi, iu

contains


!//////////////////////////////////////////////////////////////////////////////!
! File status/read/write
!//////////////////////////////////////////////////////////////////////////////!

subroutine check_error(condition,msg0,msg1,ou)
  implicit none
  logical, intent(in) :: condition
  character(len=*), intent(in) :: msg0
  character(len=*), intent(in), optional :: msg1
  integer, intent(in), optional :: ou
  integer :: u
  
  if (condition) then
    u = 6
    if (present(ou))  u = ou
    write(u,*) 'error: ', trim(msg0)
    if (present(msg1)) write(*,*), trim(msg1)
    stop
  end if

end subroutine check_error


function file_exist(fname) result(f)
  implicit none
  character(*), intent(in) :: fname
  logical :: f

  inquire(file=trim(fname),exist=f)

end function file_exist


function filecolumns_fname(f) result(n)
  implicit none
  character(*), intent(in) :: f
  integer :: n

  open(unit=20,file=trim(f),status='old')
  n = filecolumns_unit(20)
  close(20)

end function filecolumns_fname


function filecolumns_unit(u) result(n)
!* return columns in file from a given unit number
  implicit none
  integer, intent(in) :: u
  integer :: n, nh, i
  character(LEN=4096*32) :: InLine

  n  = 0
  nh = read_commentout(20)
  do i = 1, nh
    read(u,'(a)') InLine
  end do
  read(u,'(a)', end = 10) InLine
  n = TxtNumberColumns(InLine)
10 rewind u
  
end function filecolumns_unit


function filelines_fname(f) result(n)
  implicit none
  character(*), intent(in) :: f
  integer :: n

  open(unit=20,file=trim(f),status='old')
  n = filelines_unit(20)
  close(20)

end function filelines_fname


function filelines_unit(u)  result(n)
  implicit none
  integer, intent(in) :: u
  integer :: n
  character(LEN=4096) :: InLine

  n = 0
  do
    read(u,'(A)',end=10) InLine
    n = n + 1
  end do
10 rewind u
    
end function filelines_unit


function read_commentout(u)  result(f)
  implicit none
  integer, intent(in) :: u
  integer :: f
  character(128) :: strm

  f = 0
  do 
    read(u,'(A)',end=10) strm
    if (strm(1:1)=='#') f = f + 1
  end do
10 rewind(u)

end function read_commentout


function TxtNumberColumns(InLine) result(n)
!* taken from CAMB
  implicit none
  character(*) :: InLine
  integer :: n, i
  logical :: isNum    
   
  n = 0
  isNum = .false.
  do i=1, len_trim(InLIne)
    if (verify(InLine(i:i),'-+eE.0123456789') == 0) then
      if (.not. IsNum) n=n+1
      IsNum=.true.
    else
      IsNum=.false.     
    end if
  end do
   
end function TxtNumberColumns



subroutine loadtxt_1d_d(f,d1,d2,d3,d4,d5,d6,d7,fsize,rows,usecols)
  implicit none
  !I/O
  integer, intent(in), optional :: fsize(1:2), rows(:), usecols(:)
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:)
  double precision, intent(out), optional :: d2(:), d3(:), d4(:), d5(:), d6(:),d7(:)
  !internal
  integer :: n, cn, ln, r(1:2), c(1:7)
  double precision, allocatable :: dat(:,:)

  open(unit=20,file=trim(f),status='old')

  if (present(fsize)) then
    cn = fsize(1)
    ln = fsize(2)
  else
    cn = filecolumns(20)
    ln = filelines(20)
    write(*,*) 'file size is', cn, ln
  end if

  allocate(dat(cn,ln))
  read(20,*) ((dat(1:cn,n)),n=1,ln)
  close(20)

  r = [1,ln]
  if (present(rows)) r = rows
  call loadtxt_checklines(f,r,size(d1))

  c = [((n),n=1,7)]
  if (present(usecols)) c = usecols

  d1 = dat(c(1),r(1):r(2))
  if (present(d2)) d2 = dat(c(2),r(1):r(2))
  if (present(d3)) d3 = dat(c(3),r(1):r(2))
  if (present(d4)) d4 = dat(c(4),r(1):r(2))
  if (present(d5)) d5 = dat(c(5),r(1):r(2))
  if (present(d6)) d6 = dat(c(6),r(1):r(2))
  if (present(d7)) d7 = dat(c(7),r(1):r(2))

  deallocate(dat)

end subroutine loadtxt_1d_d


subroutine loadtxt_1d_c(f,dat,usecols,fsize,debug)
  implicit none
  !I/O
  logical ,intent(in), optional :: debug
  integer, intent(in), optional :: usecols(:), fsize(1:2)
  character(*), intent(in) :: f
  double complex, intent(out) :: dat(:)
  !internal
  integer :: n, cn, ln
  double precision, allocatable :: rdat(:,:)
  double complex :: r

  open(unit=20,file=trim(f),status='old')

  if (present(fsize)) then !faster if cn is given
    cn = fsize(1)
    ln = fsize(2)
  else
    cn = filecolumns(20)
    ln = filelines(20)
    write(*,*) 'file size is', cn, ln
  end if

  allocate(rdat(cn,ln))
  if (present(debug)) write(*,*) 'reading', trim(f)
  read(20,*) ((rdat(1:cn,n)),n=1,ln)
  dat = rdat(1,:) + iu*rdat(2,:)
  close(20)
  deallocate(rdat)

end subroutine loadtxt_1d_c


subroutine loadtxt_2d_d(f,d,fsize,rows,usecols,trans,debug)
  implicit none
  !I/O
  logical, intent(in), optional :: trans, debug
  integer, intent(in), optional :: fsize(2), rows(:), usecols(:)
  character(*), intent(in) :: f
  double precision, intent(out) :: d(:,:)
  !internal
  integer :: n, m, cn, ln, n1, n2, r(1:2)
  integer, allocatable :: c(:)
  double precision, allocatable :: dat(:,:)

  !* read data from file
  open(unit=20,file=trim(f),status='old')

  if (present(fsize)) then
    cn = fsize(1)
    ln = fsize(2)
  else
    cn = filecolumns(20)
    ln = filelines(20)
  end if

  if (present(debug)) write(*,*) 'reading', trim(f)
  allocate(dat(cn,ln))
  read(20,*) ((dat(1:cn,n)),n=1,ln)
  close(20)

  r = [1,ln]
  if (present(rows)) r = rows

  if (present(trans).and.trans) then
    n1 = size(d,dim=2)
    n2 = size(d,dim=1)
  else
    n1 = size(d,dim=1)
    n2 = size(d,dim=2)
  end if

  !* output
  allocate(c(n1))
  c = [((n),n=1,n1)]
  if (present(usecols)) c = usecols
  call loadtxt_checklines(f,r,n2)
  call loadtxt_checkcols(f,c,cn)
  if (present(trans).and.trans) then
    do n = 1, n1
      d(:,n) = dat(c(n),r(1):r(2))
    end do
  else
    do n = 1, n1
      d(n,:) = dat(c(n),r(1):r(2))
    end do
  end if
  deallocate(dat,c)

end subroutine loadtxt_2d_d


subroutine loadtxt_2d_c(f,d,fsize,rows,usecols,trans)
  implicit none
  !I/O
  logical, intent(in), optional :: trans
  integer, intent(in), optional :: fsize(1:2), rows(:), usecols(:)
  character(*), intent(in) :: f
  double complex, intent(out) :: d(:,:)
  !internal
  integer :: m, n, i, n1, n2, cn, ln, r(2)
  integer, allocatable :: c(:)
  double precision, allocatable :: rdat(:,:)
  double complex, allocatable :: dat(:,:)

  open(unit=20,file=trim(f),status='old')

  if (present(fsize)) then
    cn = fsize(1)
    ln = fsize(2)
  else
    cn = filecolumns(20)
    ln = filelines(20)
  end if

  allocate(dat(cn/2,ln),rdat(cn,ln))
  read(20,*) ((rdat(1:cn,n)),n=1,ln)
  do i = 1, cn/2
    dat(i,1:ln) = [((rdat(2*i-1,n)+iu*rdat(2*i,n)),n=1,ln)]
  end do
  close(20)
  deallocate(rdat)

  r = [1,ln]
  if (present(rows)) r = rows

  if (present(trans).and.trans) then
    n1 = size(d,dim=2)
    n2 = size(d,dim=1)
  else
    n1 = size(d,dim=1)
    n2 = size(d,dim=2)
  end if

  allocate(c(n1))
  c = [((n),n=1,n1)]
  if (present(usecols)) c = usecols
  call loadtxt_checklines(f,r,n2)
  call loadtxt_checkcols(f,c,cn)
  if (present(trans).and.trans) then
    do n = 1, n1
      d(:,n) = dat(c(n),r(1):r(2))
    end do
  else
    do n = 1, n1
      d(n,:) = dat(c(n),r(1):r(2))
    end do
  end if
  deallocate(dat,c)

end subroutine loadtxt_2d_c


subroutine loadtxt_checklines(f,r,n)
  implicit none
  character(*), intent(in) :: f
  integer, intent(in) :: r(1:2), n

  if (r(2)-r(1)+1/=n) then
    write(*,*) 'reading', trim(f)
    write(*,*) 'error: size of output data is inconsistent with number of file lines or specified rows'
    write(*,*) 'file lines (or specified rows):', r(2)-r(1)+1
    write(*,*) 'array size:', n
    stop
  end if

end subroutine loadtxt_checklines


subroutine loadtxt_checkcols(f,c,cn)
  implicit none
  character(*), intent(in) :: f
  integer, intent(in) :: c(:), cn

  if (maxval(c)>cn) then
    write(*,*) 'reading', trim(f)
    write(*,*) 'error: not enough file columns'
    write(*,*) 'maximum column', cn
    write(*,*) 'maximum column index', maxval(c)
    stop
  end if

end subroutine loadtxt_checkcols


subroutine ifexist_chorg(fname)
  implicit none
  character(*), intent(in) :: fname
  logical :: fexist

  inquire(file=trim(fname),exist=fexist)
  if (fexist)  call system('cp '//trim(fname)//' '//trim(fname)//'.tmp')

end subroutine ifexist_chorg


subroutine ifexist_rmorg(fname)
  implicit none
  character(*), intent(in) :: fname
  logical :: fexist

  inquire(file=trim(fname),exist=fexist)
  if (fexist)  call system('rm -rf '//trim(fname)//' '//trim(fname)//'.tmp')

end subroutine ifexist_rmorg


function fexist(fname)  result(f)
  implicit none
  character(*), intent(in) :: fname
  logical :: f

  inquire(file=trim(fname),exist=f)

end function fexist


subroutine savetxt_assign_d(d,dat,n)
!* subroutine for savetxt_1d_d
  implicit none
  integer, intent(in) :: n
  double precision, intent(in)  :: d(:)
  double precision, intent(out) :: dat(:)

  if (size(d)/=size(dat)) then
    write(*,*) 'error (savetxt): data '//str(n)//' size should be equal'
    stop
  end if
  dat = d

end subroutine savetxt_assign_d


subroutine savetxt_assign_c(d,dat,n)
!* subroutine for savetxt_1d_c
  implicit none
  integer, intent(in) :: n
  double complex, intent(in)  :: d(:)
  double complex, intent(out) :: dat(:)

  if (size(d)/=size(dat)) then
    write(*,*) 'error (savetxt): data '//str(n)//' size should be equal'
    stop
  end if
  dat = d

end subroutine savetxt_assign_c


subroutine savetxt_1d_d(f,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,ac,ow)
  implicit none
  !I/O
  logical, intent(in), optional :: ow
  character(*), intent(in) :: f
  integer, intent(in), optional :: ac
  double precision, dimension(:), intent(in) :: d1
  double precision, dimension(:), intent(in), optional :: d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16
  !internal
  character(8) :: st
  integer :: i, l, n, a0, a1, m(16)
  double precision, allocatable :: dat(:,:)

  !set output format
  a0 = 14
  if (present(ac)) a0 = ac
  a1 = a0 - 7

  ! count columns
  m    = 0
  m(1) = 1
  if (present(d2))   m(2)  = 1 + maxval(m)
  if (present(d3))   m(3)  = 1 + maxval(m)
  if (present(d4))   m(4)  = 1 + maxval(m)
  if (present(d5))   m(5)  = 1 + maxval(m)
  if (present(d6))   m(6)  = 1 + maxval(m)
  if (present(d7))   m(7)  = 1 + maxval(m)
  if (present(d8))   m(8)  = 1 + maxval(m)
  if (present(d9))   m(9)  = 1 + maxval(m)
  if (present(d10))  m(10) = 1 + maxval(m)
  if (present(d11))  m(11) = 1 + maxval(m)
  if (present(d12))  m(12) = 1 + maxval(m)
  if (present(d13))  m(13) = 1 + maxval(m)
  if (present(d14))  m(14) = 1 + maxval(m)
  if (present(d15))  m(15) = 1 + maxval(m)
  if (present(d16))  m(16) = 1 + maxval(m)

  n = maxval(m)
  l = size(d1,dim=1)

  ! set data
  allocate(dat(n,l))
  dat(m(1),:) = d1
  if (present(d2))   call savetxt_assign_d(d2,dat(m(2),:),2)
  if (present(d3))   call savetxt_assign_d(d3,dat(m(3),:),3)
  if (present(d4))   call savetxt_assign_d(d4,dat(m(4),:),4)
  if (present(d5))   call savetxt_assign_d(d5,dat(m(5),:),5)
  if (present(d6))   call savetxt_assign_d(d6,dat(m(6),:),6)
  if (present(d7))   call savetxt_assign_d(d7,dat(m(7),:),7)
  if (present(d8))   call savetxt_assign_d(d8,dat(m(8),:),8)
  if (present(d9))   call savetxt_assign_d(d9,dat(m(9),:),9)
  if (present(d10))  call savetxt_assign_d(d10,dat(m(10),:),10)
  if (present(d11))  call savetxt_assign_d(d11,dat(m(11),:),11)
  if (present(d12))  call savetxt_assign_d(d12,dat(m(12),:),12)
  if (present(d13))  call savetxt_assign_d(d13,dat(m(13),:),13)
  if (present(d14))  call savetxt_assign_d(d14,dat(m(14),:),14)
  if (present(d15))  call savetxt_assign_d(d15,dat(m(15),:),15)
  if (present(d16))  call savetxt_assign_d(d16,dat(m(16),:),16)

  ! output
  if (present(ow).and.ow==.true.) then
    st = 'replace'
  else if (present(ow).and.ow==.false.) then
    call ifexist_chorg(f)
  else
    st = 'new'
  end if
  open(unit=20,file=trim(f),status=st)
  write(20,'('//str(n)//'(E'//str(a0)//'.'//str(a1)//',1X))') ((dat(1:n,i)),i=1,l)
  close(20)

  deallocate(dat)

end subroutine savetxt_1d_d


subroutine savetxt_1d_c(f,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,ac,ow)
  implicit none
  !I/O
  logical, intent(in), optional :: ow
  character(*), intent(in) :: f
  integer, intent(in), optional :: ac
  double complex, dimension(:), intent(in) :: d1
  double complex, dimension(:), intent(in), optional :: d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16
  !internal
  character(8) :: st
  integer :: i, l, n, a0, a1, m(16)
  double complex, allocatable :: dat(:,:)

  !set output format
  a0 = 14
  if (present(ac)) a0 = ac
  a1 = a0 - 7

  ! count columns
  m    = 0
  m(1) = 1
  if (present(d2))   m(2) = 1 + maxval(m)
  if (present(d3))   m(3) = 1 + maxval(m)
  if (present(d4))   m(4) = 1 + maxval(m)
  if (present(d5))   m(5) = 1 + maxval(m)
  if (present(d6))   m(6) = 1 + maxval(m)
  if (present(d7))   m(7) = 1 + maxval(m)
  if (present(d8))   m(8)  = 1 + maxval(m)
  if (present(d9))   m(9)  = 1 + maxval(m)
  if (present(d10))  m(10) = 1 + maxval(m)
  if (present(d11))  m(11) = 1 + maxval(m)
  if (present(d12))  m(12) = 1 + maxval(m)
  if (present(d13))  m(13) = 1 + maxval(m)
  if (present(d14))  m(14) = 1 + maxval(m)
  if (present(d15))  m(15) = 1 + maxval(m)
  if (present(d16))  m(16) = 1 + maxval(m)

  n = maxval(m)
  l = size(d1,dim=1)

  ! set data
  allocate(dat(n,l))
  dat(m(1),:) = d1
  if (present(d2))   call savetxt_assign_c(d2,dat(m(2),:),2)
  if (present(d3))   call savetxt_assign_c(d3,dat(m(3),:),3)
  if (present(d4))   call savetxt_assign_c(d4,dat(m(4),:),4)
  if (present(d5))   call savetxt_assign_c(d5,dat(m(5),:),5)
  if (present(d6))   call savetxt_assign_c(d6,dat(m(6),:),6)
  if (present(d7))   call savetxt_assign_c(d7,dat(m(7),:),7)
  if (present(d8))   call savetxt_assign_c(d8,dat(m(8),:),8)
  if (present(d9))   call savetxt_assign_c(d9,dat(m(9),:),9)
  if (present(d10))  call savetxt_assign_c(d10,dat(m(10),:),10)
  if (present(d11))  call savetxt_assign_c(d11,dat(m(11),:),11)
  if (present(d12))  call savetxt_assign_c(d12,dat(m(12),:),12)
  if (present(d13))  call savetxt_assign_c(d13,dat(m(13),:),13)
  if (present(d14))  call savetxt_assign_c(d14,dat(m(14),:),14)
  if (present(d15))  call savetxt_assign_c(d15,dat(m(15),:),15)
  if (present(d16))  call savetxt_assign_c(d16,dat(m(16),:),16)

  ! output
  if (present(ow).and.ow==.true.) then
    st = 'replace'
  else if (present(ow).and.ow==.false.) then
    call ifexist_chorg(f)
  else
    st = 'new'
  end if

  open(unit=20,file=trim(f),status=st)
  write(20,'('//str(2*n)//'(E'//str(a0)//'.'//str(a1)//',1X))') ((dat(1:n,i)),i=1,l)
  close(20)

  deallocate(dat)

end subroutine savetxt_1d_c


subroutine savetxt_2d_d(f,d1,d2,d3,d4,d5,ac,ow)
  implicit none
  !I/O
  logical, intent(in), optional :: ow
  character(*), intent(in) :: f
  integer, intent(in), optional :: ac
  double precision, dimension(:,:), intent(in) :: d1
  double precision, dimension(:,:), intent(in), optional :: d2, d3, d4, d5
  !internal
  character(8) :: st
  integer :: i, l, n, m(5), a0, a1
  double precision, allocatable :: dat(:,:)

  ! set output format
  a0 = 14
  if (present(ac)) a0 = ac
  a1 = a0 - 7

  ! count columns
  m = 0
  m(1) = size(d1,dim=1)
  if (present(d2))  m(2) = size(d2,dim=1)
  if (present(d3))  m(3) = size(d3,dim=1)
  if (present(d4))  m(4) = size(d4,dim=1)
  if (present(d5))  m(5) = size(d5,dim=1)

  n = sum(m)
  l = size(d1,dim=2)

  ! set data
  allocate(dat(n,l));  dat=0d0
  dat(1:m(1),:) = d1
  if (present(d2)) then
    if (size(d2,dim=2)/=l) stop 'error: data size should be equal'
    dat(m(1)+1:sum(m(1:2)),:) = d2(1:m(2),:)
  end if
  if (present(d3)) then
    if (size(d3,dim=2)/=l) stop 'error: data size should be equal'
    dat(sum(m(1:2))+1:sum(m(1:3)),:) = d3 
  end if
  if (present(d4)) then
    if (size(d4,dim=2)/=l) stop 'error: data size should be equal'
    dat(sum(m(1:3))+1:sum(m(1:4)),:) = d4 
  end if
  if (present(d5)) then
    if (size(d5,dim=2)/=l) stop 'error: data size should be equal'
    dat(sum(m(1:4))+1:sum(m(1:5)),:) = d5 
  end if

  ! output
  if (present(ow).and.ow==.true.) then
    st = 'replace'
  else if (present(ow).and.ow==.false.) then
    call ifexist_chorg(f)
  else
    st = 'new'
  end if
  open(unit=20,file=trim(f),status=st)
  write(20,'('//str(n)//'(E'//str(a0)//'.'//str(a1)//',1X))') ((dat(1:n,i)),i=1,l)
  close(20)

  deallocate(dat)

end subroutine savetxt_2d_d


subroutine savetxt_2d_c(f,d1,d2,d3,d4,d5,ac,ow)
  implicit none
  !I/O
  logical, intent(in), optional :: ow
  character(*), intent(in) :: f
  integer, intent(in), optional :: ac
  double complex, dimension(:,:), intent(in) :: d1
  double complex, dimension(:,:), intent(in), optional :: d2, d3, d4, d5
  !internal
  character(8) :: st
  integer :: i, l, n, m(5), a0, a1
  double complex, allocatable :: dat(:,:)

  ! set output format
  a0 = 14
  if (present(ac)) a0 = ac
  a1 = a0 - 7

  ! count columns
  m = 0
  m(1) = size(d1,dim=1)
  if (present(d2))  m(2) = size(d2,dim=1)
  if (present(d3))  m(3) = size(d3,dim=1)
  if (present(d4))  m(4) = size(d4,dim=1)
  if (present(d5))  m(5) = size(d5,dim=1)

  n = sum(m)
  l = size(d1,dim=2)

  ! set data
  allocate(dat(n,l));  dat=0d0
  dat(1:m(1),:) = d1
  if (present(d2)) then
    if (size(d2,dim=2)/=l) stop 'error: data size should be equal'
    dat(m(1)+1:sum(m(1:2)),:) = d2(1:m(2),:)
  end if
  if (present(d3)) then
    if (size(d3,dim=2)/=l) stop 'error: data size should be equal'
    dat(sum(m(1:2))+1:sum(m(1:3)),:) = d3 
  end if
  if (present(d4)) then
    if (size(d4,dim=2)/=l) stop 'error: data size should be equal'
    dat(sum(m(1:3))+1:sum(m(1:4)),:) = d4 
  end if
  if (present(d5)) then
    if (size(d5,dim=2)/=l) stop 'error: data size should be equal'
    dat(sum(m(1:4))+1:sum(m(1:5)),:) = d5 
  end if

  ! output
  if (present(ow).and.ow==.true.) then
    st = 'replace'
  else if (present(ow).and.ow==.false.) then
    call ifexist_chorg(f)
  else
    st = 'new'
  end if
  open(unit=20,file=trim(f),status=st)
  write(20,'('//str(2*n)//'(E'//str(a0)//'.'//str(a1)//',1X))') ((dat(1:n,i)),i=1,l)
  close(20)

  deallocate(dat)

end subroutine savetxt_2d_c


!//////////////////////////////////////////////////////////////////////////////!
! interface linspace 
!//////////////////////////////////////////////////////////////////////////////!

function linspace_dble(ini,fin,n)  result(f)
!* output ini, ini*d, ini*2*d, ..., fin
  implicit none
!
! [input]
!   ini --- initial value
!   fin --- final value
!   n   --- number of points
  integer, intent(in) :: n
  double precision, intent(in) :: ini, fin
!
! [internal]
  integer :: i
  double precision :: f(1:n), d

  d = (fin-ini)/dble(n-1)
  f = [ ( (ini+d*(i-1)), i=1,n ) ]

end function linspace_dble


function linspace_dble_div(div,n)  result(f)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: div(1:2)
  double precision :: f(1:n)

  f = linspace_dble(div(1),div(2),n)

end function linspace_dble_div


function linspace_int(ini,fin,n)  result(f)
  implicit none
  integer, intent(in) :: n, ini, fin
  integer :: i
  double precision :: f(1:n), d

  d = dble(fin-ini)/dble(n-1)
  f = [ ( (ini+d*(i-1)), i=1,n ) ]

end function linspace_int


function linspace_int_simple(ini,n)  result(f)
  implicit none
  integer, intent(in) :: ini, n
  integer :: i
  double precision :: f(1:n-ini+1)

  f = [ ( (ini+i-1), i=1,n-ini+1 ) ]

end function linspace_int_simple


function linspace_int_div(div,n)  result(f)
  implicit none
  integer, intent(in) :: n, div(1:2)
  double precision :: f(1:n)

  f = linspace_int(div(1),div(2),n)

end function linspace_int_div


!//////////////////////////////////////////////////////////////////////////////!
! interpolation 
!//////////////////////////////////////////////////////////////////////////////!

function interp_lin_narr(x,x0,x1,y0,y1)  result(f)
  implicit none
  double precision, intent(in) :: x,x0,x1,y0,y1
  double precision :: f

  f = y0 + ((x-x0)/(x1-x0))*(y1-y0)

end function interp_lin_narr


function interp_lin_arr(x,xn,yn) result(f)
  implicit none
  double precision, intent(in) :: x(:),xn(:),yn(:)
  double precision :: f(size(x))
  integer :: i, n
  double precision :: x0, x1, y0, y1

  do i = 1, size(x)
    n    = neighb(x(i),xn)
    x0   = xn(n)
    x1   = xn(n+1)
    y0   = yn(n)
    y1   = yn(n+1)
    f(i) = y0 + ((x(i)-x0)/(x1-x0))*(y1-y0)
  end do

end function interp_lin_arr


subroutine spline(x,y,n,yp1,ypn,y2)
! * from numerical recipe, modified for f90.
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n. 
  implicit none
  !I/O
  integer, intent(in) :: n
  double precision, intent(in)  :: yp1, ypn, x(n), y(n)
  double precision, intent(out) :: y2(n)
  !internal
  integer :: i, k
  double precision :: p, qn, sig, un, u(n)

  if (yp1.gt..99d30) then
    y2(1) = 0d0
    u(1)  = 0d0
  else
    y2(1) = -0.5d0
    u(1)  = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  end if

  do i=2, n-1
    sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p     = sig*y2(i-1)+2d0
    y2(i) = (sig-1d0)/p
    u(i)  = (6d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) / (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do

  if (ypn.gt..99d30) then
    qn = 0d0
    un = 0d0
  else
    qn = 0.5d0
    un = (3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1d0)
  do k = n-1, 1, -1
    y2(k) = y2(k)*y2(k+1) + u(k)
  end do

end subroutine spline 


function splint(x,xa,ya,y2a)  result(f)
  implicit none
  double precision, intent(in) :: x,xa(:),y2a(:),ya(:)
  double precision :: f
  integer :: khi,klo

  klo = neighb(x,xa)
  khi = klo + 1
  f = interp_sp(x,xa(klo),xa(khi),ya(klo),ya(khi),y2a(klo),y2a(khi))

end function splint


function neighb(x,arr)  result(nn)
! return nearest neighbour (nn) which satisfies arr(nn)<=x<=arr(nn+1)
  implicit none
  double precision, intent(in) :: x, arr(:)
  integer :: nn, n, lo, hi

  lo = 1
  hi = size(arr)

1 if (hi-lo.gt.1) then
    n = (hi+lo)/2
    if(arr(n).gt.x)then
      hi = n
    else
      lo = n
    end if
    goto 1
  end if
  nn = lo

end function neighb


function interp_sp(x,x0,x1,y0,y1,sp0,sp1)  result(f)
  implicit none
  double precision, intent(in) :: x,x0,x1,y0,y1,sp0,sp1
  double precision :: f, a, b, h

  h = x1-x0
  if (h.eq.0.) stop 'error (interp_sp): h is zero'
  a = (x1-x)/h
  b = (x-x0)/h
  f = a*y0 + b*y1 + ((a**3-a)*sp0+(b**3-b)*sp1)*h**2/6.0
 
end function interp_sp


!//////////////////////////////////////////////////////////////////////////////!
! MC simulation 
!//////////////////////////////////////////////////////////////////////////////!

subroutine meanvar(Cl,mCl,vCl,f) 
!Compute mean and variance of 1D array
!Input: 1D array for number of realizations (Cl)
!Outputs: mean value (mCl), variance (vCl) 
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  double precision, intent(in) :: Cl(:,:)
  double precision, intent(out), optional :: mCl(:), vCl(:)
  !internal
  integer :: n, l, i
  double precision, allocatable :: mean(:), var(:)

  n = size(Cl,dim=1)
  l = size(Cl,dim=2)

  allocate(mean(l),var(l)); mean=0d0; var=0d0

  if(n == 1) then 
    mean = Cl(1,:)
  else if (n >= 2) then 
    mean = sum(Cl,dim=1)/dble(n) !mean
    var  = dsqrt( sum(Cl**2,dim=1)/dble(n) - mean**2 ) !variance
  end if

  if(present(mCl)) mCl = mean
  if(present(vCl)) vCl = var
  if(present(f))   call savetxt(f,linspace(1,l,l),mean,var)
  deallocate(mean,var)

end subroutine meanvar


function ave_1d(x)  result(f)
  implicit none
  double precision, intent(in) :: x(:)
  double precision, allocatable :: f

  f = sum(x)/dble(size(x))

end function ave_1d


function ave_2d(x,d)  result(f)
  implicit none
  integer, intent(in) :: d
  double precision, intent(in) :: x(:,:)
  double precision, allocatable :: f(:)
  integer :: n, s

  if (d==1) n = size(x,dim=2)
  if (d==2) n = size(x,dim=1)

  s = size(x,dim=d)

  allocate(f(n))
  f = sum(x,dim=d)/dble(s)


end function ave_2d


subroutine save_average_2d(f,dat,id,bc)
!average 2d quantities along an axis (1st index for default), compute its mean and variance, and save to a file
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(in) :: dat(:,:,:)
  !(optional)
  integer, intent(in), optional :: id(2)
  double precision, intent(in), optional :: bc(:)
  !internal
  integer :: ndat, bmax, j
  double precision, allocatable :: b(:,:), mdat(:,:), vdat(:,:)

  if (present(id)) then
    ndat = size(dat,id(1))
    bmax = size(dat,id(2))
  else
    ndat = size(dat,2)
    bmax = size(dat,3)
  end if

  if (ndat<1) stop 'save_average_2d: ndat is strange'
  if (bmax<1) stop 'save_average_2d: bmax is strange'

  allocate(mdat(ndat,bmax),vdat(ndat,bmax),b(1,bmax)); mdat=0d0; vdat=0d0
  do j = 1, ndat
    call meanvar(dat(:,j,:),mdat(j,:),vdat(j,:))
  end do

  write(*,*) 'output average and variance'
  if (present(bc)) then
    b(1,:) = bc
    call savetxt(f,b,mdat,vdat,ow=.true.)
  else
    b(1,:) = linspace(1,bmax)
    call savetxt(f,mdat,vdat,ow=.true.)
  end if

  deallocate(mdat,vdat)

end subroutine save_average_2d


!//////////////////////////////////////////////////////////////////////////////!
! histogram 
!//////////////////////////////////////////////////////////////////////////////!

subroutine get_histogram_bin(a1,a2,bn,bp)
  implicit none
  !I/O
  integer, intent(in) :: bn
  double precision, intent(in) :: a1, a2
  double precision, intent(out) :: bp(0:bn)
  !internal
  integer :: b
  double precision :: w

  w = (a2-a1)/dble(bn)
  do b = 0, bn
    bp(b) = dble(b)*w + a1
  end do

end subroutine get_histogram_bin


subroutine get_histogram(hist,a1,a2,f,x)
  implicit none
  !I/O
  double precision, intent(in) :: a1, a2, f(:)
  double precision, intent(in), optional :: x(:)
  double precision, intent(out) :: hist(:)
  !internal
  integer :: b, i, xn, bn
  double precision, allocatable :: bp(:), k(:), h(:)

  bn = size(hist)
  allocate(bp(0:bn))

  if (present(x)) then
    xn = size(x)
    allocate(k(xn),h(xn))
    k  = x
    h  = f
  else
    xn = size(f)
    allocate(k(xn),h(xn))
    k  = f
    h  = 1d0
  end if

  call get_histogram_bin(a1,a2,bn,bp)

  !make histogram
  hist = 0
  do b = 1, bn
    do i = 1, xn
      if(bp(b-1)<k(i).and.k(i)<=bp(b)) hist(b) = hist(b) + h(i)
    end do
  end do
  hist = hist/dble(xn)

  deallocate(bp,k,h)


end subroutine get_histogram

!//////////////////////////////////////////////////////////////////////////////!
! Quadratures 
!//////////////////////////////////////////////////////////////////////////////!

subroutine get_zeropoint(GLn,GLw,GLz,GLeps)
! Zero points of Gauss-Legendre !
  implicit none
  !I/O
  integer, intent(in) :: GLn
  double precision, intent(in) :: GLeps
  double precision, intent(out) :: GLw(:), GLz(:)
  !internal
  integer :: m, j, i
  double precision :: z1, z, pp, p3, p2 ,p1

  m = (GLn+1)/2
  do i = 1, m
    z = dcos(pi*(dble(i)-0.25d0)/(dble(GLn)+0.5d0))
1   continue
    p1 = 1d0
    p2 = 0d0
    do j = 1, GLn
      p3 = p2
      p2 = p1
      p1 = ((2d0*j-1d0)*z*p2-(j-1d0)*p3)/dble(j)
    end do
    pp = GLn*(z*p1-p2)/(z*z-1d0)
    z1 = z
    z  = z1-p1/pp
    if(abs(z-z1).gt.GLeps) goto 1
    GLz(i) = -z
    GLz(GLn+1-i) = z
    GLw(i) = 2d0/((1d0-z*z)*pp*pp)
    GLw(GLn+1-i) = GLw(i)
  end do

end subroutine get_zeropoint


subroutine gl_initialize(GL,gln,eps)
  implicit none
  integer, intent(in) :: gln
  double precision, intent(in), optional :: eps
  type(GAUSS_LEGENDRE_PARAMS), intent(out) :: GL

  GL%n = gln
  GL%EPS = 5d-16
  if (present(eps)) GL%EPS = eps
  allocate(GL%z(GL%n),GL%w(GL%n))
  call get_zeropoint(GL%n,GL%w,GL%z,GL%eps)

end subroutine gl_initialize


subroutine gl_finalize(GL)
  implicit none
  type(GAUSS_LEGENDRE_PARAMS), intent(inout) :: GL

  deallocate(GL%z,GL%w)

end subroutine gl_finalize


subroutine ZeroPoint(GL)
  implicit none
  !I/O
  type(GAUSS_LEGENDRE_PARAMS), intent(inout) :: GL

  call get_zeropoint(GL%n,GL%w,GL%z,GL%eps)

end subroutine ZeroPoint


function GLpoint(xran,GLz)  result(x)
!* x values for GL quadrature
!    (b+a)/2 + (b-a)/2 x_{i,0}
! where x_{i,0} are the solution of P_{n+1}(x)=0
  implicit none
  double precision :: xran(:), GLz
  double precision :: x, xm, xl

  xm = 0.5d0*(xran(2)+xran(1))
  xl = 0.5d0*(xran(2)-xran(1))
  x = xm + xl*GLz

end function GLpoint


function GLpoints(xran,GLz)  result(x)
  implicit none
  double precision :: xran(:), GLz(:)
  double precision :: x(size(GLz))
  integer :: i

  do i = 1, size(GLz)
    x(i) = GLpoint(xran,GLz(i))
  end do

end function GLpoints


function GLdx(xran,GLw)  result(dx)
!* interval for GL quadrature
!    (b-a)/2 * w_i
  implicit none
  double precision :: xran(:), GLw
  double precision :: dx, xl

  xl = 0.5d0*(xran(2)-xran(1))
  dx = xl*GLw

end function GLdx


function GLdxs(xran,GLw)  result(dx)
  implicit none
  double precision :: xran(:), GLw(:)
  double precision :: dx(size(GLw)), xl
  integer :: i

  do i = 1, size(GLw)
    dx(i) = GLdx(xran,GLw(i))
  end do

end function GLdxs


!///////////////////////////////////////////////////////!
! * Matrix Analysis
!///////////////////////////////////////////////////////!

subroutine multigrid_obscov(Cov,CovM,npix,mpix,s,nn,symmetric)
!* transform pixel-pixel covariance into multigride space
  implicit none
  logical, intent(in), optional :: symmetric
  integer, intent(in) :: npix, mpix, nn(1:2), s
  double precision, intent(in) :: Cov(npix,npix)
  double precision, intent(out) :: CovM(mpix,mpix)
  !internal
  logical :: sym = .false.
  integer, allocatable :: label(:,:,:), a(:,:), b(:,:)
  integer :: mi, mj, si, sj
  double precision :: x

  if(present(symmetric)) sym = symmetric
  allocate(label(0:mpix-1,s,2))
  call multigrid_index(nn(1),mpix,s,label)

  allocate(a(s,2),b(s,2))
  do mi = 1, mpix
    do mj = 1, mpix
      if(sym.and.mj<mi) cycle
      a = label(mi-1,:,:)
      b = label(mj-1,:,:)
      if(mi==1.and.mj==1) write(*,*) a(:,:)
      x = 0d0
      do si = 1, s
        do sj = 1, s
          if(product(Cov(a(si,1):a(si,2),b(sj,1):b(sj,2)))==0d0) then
            x = 0d0
            goto 11
          end if
          x = x + sum(Cov(a(si,1):a(si,2),b(sj,1):b(sj,2)))
        end do
      end do
11    CovM(mi,mj) = x/dble(s**4)
      if(sym.and..not.mi==mj) CovM(mj,mi) = CovM(mi,mj)
    end do
  end do  
  deallocate(label,a,b)


end subroutine multigrid_obscov


subroutine multigrid_index(mmx,mpix,s,label)
!generate matrix index for multigrid
  implicit none
  !I/O
  integer, intent(in) :: mmx,mpix,s
  integer, intent(out) :: label(0:mpix-1,s,2)
  !internal
  integer :: m, xi, yi, si, pix

  do m = 0, mpix-1
    do si = 1, s
      xi = mod(m,mmx)
      yi = int(m/mmx)
      pix = yi*s**2*mmx
      label(m,si,1) = pix + 1 + s*xi + (si-1)*s*mmx
      label(m,si,2) = pix + s*(xi+1) + (si-1)*s*mmx
    end do
  end do

end subroutine multigrid_index


function vec2mat_dble(a,b) result(M)
  implicit none
  double precision, intent(in) :: a(:), b(:)
  double precision :: M(size(a),size(b))
  integer :: i,j

  do i = 1, size(a)
    do j = 1, size(b)
      M(i,j) = a(i)*b(j)
    end do
  end do

end function vec2mat_dble


function vec2mat_cmplx(a,b) result(M)
  implicit none
  double complex, intent(in) :: a(:), b(:)
  double complex :: M(size(a),size(b))
  integer :: i,j

  do i = 1, size(a)
    do j =1, size(b)
      M(i,j) = a(i)*b(j)
    end do
  end do

end function vec2mat_cmplx


function vec2mat_diag(vec) result(M)
  implicit none
  double precision, intent(in) :: vec(:)
  integer :: i
  double precision :: M(size(vec),size(vec))

  M = 0d0
  do i = 1, size(vec)
    M(i,i) = vec(i)
  end do

end function vec2mat_diag


function mdiag(i,n)
!* diagonal elements of n*n matrix
! 1 + n + (n-1) + (n-2) + ...
!   = n*(n+1)/2 - (n-i)(n-i+1)/2 - (n-i)
  implicit none
  !I/O
  integer, intent(in) :: i, n
  integer :: MDIAG

  MDIAG = n*(i-1) + i*(3-i)/2

end function mdiag


function mlab(i,j,n)
!* elements of n*n matrix
!(1,1), (1,2), (1,3), ..., (2,2), (2,3), ... -> 1,2,3,...
  implicit none
  !I/O
  integer, intent(in) :: i,j,n
  integer :: MLAB

  if(j>=i) MLAB = mdiag(i,n) + (j-i)
  if(j<i)  MLAB = mdiag(j,n) + (i-j)

end function mlab


subroutine symmetric(M)
! symmetrize matrix (M_ji -> M_ij)
  implicit none
  !I/O
  double precision, intent(inout) :: M(:,:)
  !internal
  integer :: i, j, n

  n = size(M,dim=1)
  do i = 1, n
    do j = i + 1, n
      M(j,i) = M(i,j)
    end do
  end do

end subroutine symmetric


function mat_identity(N) result(M)
! get N*N indentity matrix
  implicit none
  integer, intent(in) :: N
  integer :: i
  double precision :: M(N,N)

  M = 0d0
  do i = 1, N
    M(i,i) = 1d0
  end do

end function mat_identity


subroutine matrix_diag_tri(n,cl,A)
  !obtain a matrix A satisfying AA^t = cl
  implicit none
  !I/O
  integer, intent(in) :: n
  double precision, intent(in), dimension(n,n) :: cl
  double precision, intent(out), dimension(n,n) :: A
  !internal
  integer :: i, k

  A = 0d0
  do k = 1, n
    ! A(k,i) = 0 if k<i
    do i = 1, k -1
      A(k,i) = ( cl(i,k) - sum(A(i,1:i-1)*A(k,1:i-1)) ) / A(i,i)
    end do
    A(k,k) = dsqrt( cl(k,k) - sum(A(k,1:k-1)*A(k,1:k-1)) )
  end do

end subroutine matrix_diag_tri



!//////////////////////////////////////////////////////////////////////////////!
! Re-ordering array
!//////////////////////////////////////////////////////////////////////////////!

subroutine sort_1d(refval,ii)
!Sorting input array 
  implicit none
  !I/O
  integer, intent(out), optional :: ii(:) !label
  double precision, intent(inout) :: refval(:) !values
  !internal
  integer :: i, n, m, id, pn
  double precision :: dummy
  double precision, allocatable :: array(:)

  pn = size(refval)
  allocate(array(pn))

  !set array
  array = refval
  if (present(ii)) then
    do n = 1, pn
      ii(n) = n
    end do
  end if

  !sort
  do n = 1, pn
    do m = n + 1, pn
      if(array(n)>=array(m)) cycle
      dummy = array(n)
      array(n) = array(m)
      array(m) = dummy
      if (present(ii)) then
        id = ii(n)
        ii(n) = ii(m)
        ii(m) = id
      end if
    end do
  end do 

  refval = array

  deallocate(array)

end subroutine sort_1d


subroutine sort_2d(refval,ii,jj)
  implicit none
  !I/O
  integer, intent(inout) :: ii(:), jj(:)
  double precision, intent(inout) :: refval(:,:)
  !internal
  integer :: i, j, n, m, id, jd, pn1, pn2
  double precision :: dummy
  double precision, allocatable :: array(:)

  pn1 = size(refval,dim=1)
  pn2 = size(refval,dim=2)

  allocate(array(pn1*pn2))

  !set to array
  n = 1
  do i = 1, pn1
    do j = 1, pn2
      array(n) = refval(i,j) 
      ii(n) = i
      jj(n) = j
      n = n + 1
    end do
  end do

  !sort
  do n = 1, pn1*pn2
    do m = n + 1, pn1*pn2
      if(array(n)>=array(m)) cycle
      dummy = array(n)
      id = ii(n)
      jd = jj(n)
      array(n) = array(m)
      ii(n) = ii(m)
      jj(n) = jj(m)
      array(m) = dummy
      ii(m) = id
      jj(m) = jd
    end do
  end do

  do n=1, pn1*pn2
    refval(ii(n),jj(n)) = array(n)
  end do

  deallocate(array)

end subroutine sort_2d



!//////////////////////////////////////////////////////////////////////////////!
! other utils
!//////////////////////////////////////////////////////////////////////////////!

subroutine check_positive(dat,p)
! check whether dat is all positive
  implicit none
  double precision, intent(in) :: dat(:)
  logical, intent(out) :: p
  integer :: i

  p = .true.
  do i = 1, size(dat)
    if (dat(i)<0d0) p = .false.
  end do

end subroutine check_positive


function param_id(str,LAB)
  implicit none
  character(*), intent(in) :: str, LAB(:)
  integer :: param_id
  integer :: p

  do p = 1, size(LAB)
    if(trim(LAB(p))==str) param_id = p
  end do

end function param_id


function str_n(n)  result(f)
  implicit none
  integer, intent(in) :: n
  character(1+idint(log10(dble(n+0.9)))) :: f
  character(LEN=128) :: a

  write(a,*) n
  f = adjustl(a)

end function str_n


function str_fix(n,l)  result(f)
  implicit none 
  integer, intent(in) :: n, l
  character(l) :: f

  write(f,'(i'//str(l)//'.'//str(l)//')') n

end function str_fix


function str_fix_r(n,l)  result(f)
  implicit none
  real, intent(in) :: n
  integer, intent(in) :: l
  character(l) :: f
  integer :: m

  m = floor(log10(n))

  if (l<abs(m)+2) stop 'error in str_fix_r: not enough length specified'
  if (n<1)  write(f,'(f'//str_n(l)//'.'//str_n(l-2)//')')   n
  if (n>=1) write(f,'(f'//str_n(l)//'.'//str_n(l-2-m)//')') n

end function str_fix_r


function str_fix_d(n,l)  result(f)
  implicit none
  double precision, intent(in) :: n
  integer, intent(in) :: l
  character(l) :: f
  integer :: m

  m = floor(log10(n))

  if (l<abs(m)+2) stop 'error in str_fix_d: not enough length specified'
  if (n<1)  write(f,'(f'//str_n(l)//'.'//str_n(l-2)//')')   n
  if (n>=1) write(f,'(f'//str_n(l)//'.'//str_n(l-2-m)//')') n

end function str_fix_d


function ones_i(n)  result(f)
  implicit none
  integer, intent(in) :: n
  integer :: i, f(n)

  f = [ ( (1), i=1,n ) ]

end function ones_i


function ones_d(n)  result(f)
  implicit none
  integer, intent(in) :: n
  integer :: i
  double precision :: f(n)

  f = [ ( (1d0), i=1,n ) ]

end function ones_d


function zeros_i(n)  result(f)
  implicit none
  integer, intent(in) :: n
  integer :: i, f(n)

  f = [ ( (0), i=1,n ) ]

end function zeros_i


function zeros_d(n)  result(f)
  implicit none
  integer, intent(in) :: n
  integer :: i
  double precision :: f(n)

  f = [ ( (0d0), i=1,n ) ]

end function zeros_d


end module general


