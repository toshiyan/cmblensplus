!////////////////////////////////////////////////////!
! * I/O 
!////////////////////////////////////////////////////!

module io
  implicit none
  integer, parameter :: dlc = KIND(0d0)
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0
  complex(dlc), parameter :: iu = (0d0,1d0)

  interface loadfile
    module procedure loadfile_1d_d, loadfile_1d_c, loadfile_2d_d, loadfile_2d_c, loadfile_3d_d, loadfile_3d_c, loadfile_4d_c
  end interface loadfile

  interface loadfile_py
    module procedure loadfile_1d_d_py, loadfile_2d_d_py
  end interface loadfile_py

  interface savefile
    module procedure savefile_1d_d, savefile_1d_c, savefile_2d_d, savefile_2d_c, savefile_3d_d, savefile_3d_c, savefile_4d_c
  end interface savefile

  private dlc, pi, iu

contains

subroutine loadfile_1d_d(f,d1,d2,d3,d4,d5,d6,d7,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:)
  double precision, intent(out), optional :: d2(:), d3(:), d4(:), d5(:), d6(:),d7(:)
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac
  integer :: frm

  ac = 'stream'
  if (present(a)) ac=a

  open(unit=20,file=trim(f),status='old',form='unformatted',access=ac)
  read(20) frm !get size

  if (size(d1)>frm) then
    write(*,*) 'error (loadfile_1d_d): size is strange'
    write(*,*) 'size specified:', size(d1)
    write(*,*) 'size read from file:', frm
    stop
  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7
  close(20)

end subroutine loadfile_1d_d


subroutine loadfile_1d_c(f,d1,d2,d3,d4,d5,d6,d7,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:)
  complex(dlc), intent(out), optional :: d2(:), d3(:), d4(:), d5(:), d6(:),d7(:)
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac
  integer :: frm

  ac = 'stream'
  if (present(a)) ac=a

  open(unit=11,file=trim(f),status='old',form='unformatted',access=ac)
  read(11) frm !get size

  if (size(d1)>frm) then
    write(*,*) 'error (loadfile_1d_c): size is strange'
    write(*,*) 'size specified:', size(d1)
    write(*,*) 'size read from file:', frm
    stop
  end if

  read(11) d1
  if (present(d2)) read(11) d2
  if (present(d3)) read(11) d3
  if (present(d4)) read(11) d4
  if (present(d5)) read(11) d5
  if (present(d6)) read(11) d6
  if (present(d7)) read(11) d7
  close(11)

end subroutine loadfile_1d_c


subroutine loadfile_1d_d_py(f,d1)
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:)
  !internal
  integer :: dummy

  open(unit=20,file=trim(f),status='old',form='unformatted',access='stream')
  read(20) dummy, d1
  close(20)

end subroutine loadfile_1d_d_py


subroutine loadfile_2d_d(f,d1,d2,d3,d4,d5,d6,d7,ns,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:,:)
  double precision, intent(out), dimension(:,:), optional :: d2,d3,d4,d5,d6,d7
  character(*), intent(in), optional :: a
  logical, intent(in), optional :: ns
  !internal
  character(16) :: ac
  integer :: i, frm(2)

  ac = 'stream'
  if (present(a)) ac=a

  !* read data from file
  open(unit=20,file=trim(f),status='old',form='unformatted',access=ac)

  if (.not.present(ns).or.ns==.false.) then

    read(20) frm !get row and column

    do i = 1, 2
      if (size(d1,dim=i)>frm(i)) then
        write(*,*) 'error (loadfile_2d_d): size is strange'
        write(*,*) 'dim:', i
        write(*,*) 'size specified:', size(d1,dim=i)
        write(*,*) 'size read from file:', frm(i)
        stop
      end if
    end do

  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7

  close(20)

end subroutine loadfile_2d_d


subroutine loadfile_2d_c(f,d1,d2,d3,d4,d5,d6,d7,ns,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:,:)
  complex(dlc), intent(out), dimension(:,:), optional :: d2,d3,d4,d5,d6,d7
  !optional
  character(*), intent(in), optional :: a
  logical, intent(in), optional :: ns
  !internal
  character(16) :: ac
  integer :: i, frm(1:2)

  ac = 'stream'
  if (present(a)) ac=a

  !* read data from file
  open(unit=20,file=trim(f),status='old',form='unformatted',access=ac)

  if (.not.present(ns).or.ns==.false.) then

    read(20) frm(1), frm(2) !get row and column

    do i = 1, 2
      if (size(d1,dim=i)>frm(i)) then
        write(*,*) 'error (loadfile_2d_c): size is strange'
        write(*,*) 'dim:', i
        write(*,*) 'size specified:', size(d1,dim=i)
        write(*,*) 'size read from file:', frm(i)
        stop
      end if
    end do

  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7

  close(20)

end subroutine loadfile_2d_c


subroutine loadfile_2d_d_py(f,d1,d2,d3,d4,d5,d6,d7)
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:,:)
  double precision, intent(out), dimension(:,:), optional :: d2,d3,d4,d5,d6,d7
  !internal
  character(16) :: ac
  integer :: i, dummy

  !* read data from file
  open(unit=20,file=trim(f),status='old',form='unformatted',access='stream')
  read(20) dummy, d1
  if (present(d2)) read(20) dummy, dummy, d2
  if (present(d3)) read(20) dummy, dummy, d3
  if (present(d4)) read(20) dummy, dummy, d4
  if (present(d5)) read(20) dummy, dummy, d5
  if (present(d6)) read(20) dummy, dummy, d6
  if (present(d7)) read(20) dummy, dummy, d7
  close(20)

end subroutine loadfile_2d_d_py


subroutine loadfile_3d_d(f,d1,d2,d3,d4,d5,d6,d7,ns,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:,:,:)
  double precision, intent(out), dimension(:,:,:), optional :: d2,d3,d4,d5,d6,d7
  !optional
  character(*), intent(in), optional :: a
  logical, intent(in), optional :: ns
  !internal
  character(16) :: ac
  integer :: i, frm(1:3)

  ac = 'stream'
  if (present(a)) ac=a

  !* read data from file
  open(unit=20,file=trim(f),status='old',form='unformatted',access=ac)

  if (.not.present(ns).or.ns==.false.) then
    read(20) frm !size
    call check_inputdata_size(frm,shape(d1),'loadfile_3d_d')
  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7

  close(20)

end subroutine loadfile_3d_d


subroutine loadfile_3d_c(f,d1,d2,d3,d4,d5,d6,d7,ns,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:,:,:)
  complex(dlc), intent(out), dimension(:,:,:), optional :: d2,d3,d4,d5,d6,d7
  !optional
  character(*), intent(in), optional :: a
  logical, intent(in), optional :: ns
  !internal
  character(16) :: ac
  integer :: i, frm(1:3)

  ac = 'stream'
  if (present(a)) ac=a

  !* read data from file
  open(unit=20,file=trim(f),status='old',form='unformatted',access=ac)

  if (.not.present(ns).or.ns==.false.) then
    read(20) frm(1), frm(2), frm(3) !get size
    call check_inputdata_size(frm,shape(d1),'loadfile_3d_c')
  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7

  close(20)

end subroutine loadfile_3d_c


subroutine loadfile_4d_c(f,d1,d2,d3,d4,d5,ns,a)
  implicit none
  !I/O
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:,:,:,:)
  complex(dlc), intent(out), dimension(:,:,:,:), optional :: d2,d3,d4,d5
  !optional
  character(*), intent(in), optional :: a
  logical, intent(in), optional :: ns
  !internal
  character(16) :: ac
  integer :: i, frm(1:3)

  ac = 'stream'
  if (present(a)) ac=a

  !* read data from file
  open(unit=20,file=trim(f),status='old',form='unformatted',access=ac)

  if (.not.present(ns).or.ns==.false.) then
    read(20) frm(1), frm(2), frm(3), frm(4) !get size
    call check_inputdata_size(frm,shape(d1),'loadfile_4d_c')
  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5

  close(20)

end subroutine loadfile_4d_c


subroutine check_inputdata_size(frm,dsize,text)
  implicit none
  character(*), intent(in) :: text
  integer, intent(in) :: frm(:), dsize(:)
  integer :: i

  do i = 1, size(frm)
    if (dsize(i)>frm(i)) then
      write(*,*) 'error ('//trim(text)//'): size is strange'
      write(*,*) 'dim:', i
      write(*,*) 'size you specified:', dsize(i)
      write(*,*) 'size read from file:', frm(i)
      stop
    end if
  end do

end subroutine check_inputdata_size


subroutine savefile_1d_d(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  double precision, dimension(:), intent(in) :: d1
  double precision, dimension(:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_1d_d


subroutine savefile_1d_c(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), dimension(:), intent(in) :: d1
  complex(dlc), dimension(:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_1d_c


subroutine savefile_2d_d(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  double precision, dimension(:,:), intent(in) :: d1
  double precision, dimension(:,:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1,dim=1), size(d1,dim=2)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_2d_d


subroutine savefile_2d_c(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), dimension(:,:), intent(in) :: d1
  complex(dlc), dimension(:,:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1,dim=1), size(d1,dim=2)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_2d_c


subroutine savefile_3d_d(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  double precision, dimension(:,:,:), intent(in) :: d1
  double precision, dimension(:,:,:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1,dim=1), size(d1,dim=2), size(d1,dim=3)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_3d_d


subroutine savefile_3d_c(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), dimension(:,:,:), intent(in) :: d1
  complex(dlc), dimension(:,:,:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1,dim=1), size(d1,dim=2), size(d1,dim=3)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_3d_c


subroutine savefile_4d_c(f,d1,d2,d3,d4,d5,ow,ns,a)
  implicit none
  !I/O
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), dimension(:,:,:,:), intent(in) :: d1
  complex(dlc), dimension(:,:,:,:), intent(in), optional :: d2, d3, d4, d5
  character(*), intent(in), optional :: a
  !internal
  character(16) :: ac, st

  !access
  ac = 'stream'
  if (present(a)) ac=a

  ! ow
  st = 'new'
  if (present(ow).and.ow==.true.) st = 'replace'

  ! output
  open(unit=20,file=trim(f),status=st,form='unformatted',access=ac)
  if (.not.present(ns).or.ns==.false.)  write(20) size(d1,dim=1), size(d1,dim=2), size(d1,dim=3), size(d1,dim=4)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)

end subroutine savefile_4d_c


end module io

