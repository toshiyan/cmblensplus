!////////////////////////////////////////////////////!
! * I/O
!////////////////////////////////////////////////////!

module io
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dlc = real64

  interface loadfile
    module procedure loadfile_1d_d, loadfile_1d_c, &
                     loadfile_2d_d, loadfile_2d_c, &
                     loadfile_3d_d, loadfile_3d_c, &
                     loadfile_4d_c
  end interface loadfile

  interface loadfile_py
    module procedure loadfile_1d_d_py, loadfile_2d_d_py
  end interface loadfile_py

  interface savefile
    module procedure savefile_1d_d, savefile_1d_c, &
                     savefile_2d_d, savefile_2d_c, &
                     savefile_3d_d, savefile_3d_c, &
                     savefile_4d_c
  end interface savefile

  private :: dlc

contains

function access_mode(a) result(ac)
  implicit none
  character(*), intent(in), optional :: a
  character(16) :: ac

  ac = 'stream'
  if (present(a)) ac = a
end function access_mode

function write_size_header(ns) result(write_size)
  implicit none
  logical, intent(in), optional :: ns
  logical :: write_size

  write_size = .true.
  if (present(ns)) write_size = .not. ns
end function write_size_header

function file_status(ow) result(st)
  implicit none
  logical, intent(in), optional :: ow
  character(16) :: st

  st = 'new'
  if (present(ow)) then
    if (ow) st = 'replace'
  end if
end function file_status

subroutine check_inputdata_size(frm, dsize, text)
  implicit none
  character(*), intent(in) :: text
  integer, intent(in) :: frm(:), dsize(:)
  integer :: i

  do i = 1, size(frm)
    if (dsize(i) > frm(i)) then
      write(*,*) 'error ('//trim(text)//'): size is strange'
      write(*,*) 'dim:', i
      write(*,*) 'size you specified:', dsize(i)
      write(*,*) 'size read from file:', frm(i)
      stop
    end if
  end do
end subroutine check_inputdata_size

subroutine loadfile_1d_d(f, d1, d2, d3, d4, d5, d6, d7, a)
  implicit none
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:)
  double precision, intent(out), optional :: d2(:), d3(:), d4(:), d5(:), d6(:), d7(:)
  character(*), intent(in), optional :: a
  integer :: frm

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  read(20) frm
  call check_inputdata_size([frm], [size(d1)], 'loadfile_1d_d')

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7
  close(20)
end subroutine loadfile_1d_d

subroutine loadfile_1d_c(f, d1, d2, d3, d4, d5, d6, d7, a)
  implicit none
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:)
  complex(dlc), intent(out), optional :: d2(:), d3(:), d4(:), d5(:), d6(:), d7(:)
  character(*), intent(in), optional :: a
  integer :: frm

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  read(20) frm
  call check_inputdata_size([frm], [size(d1)], 'loadfile_1d_c')

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  if (present(d6)) read(20) d6
  if (present(d7)) read(20) d7
  close(20)
end subroutine loadfile_1d_c

subroutine loadfile_1d_d_py(f, d1)
  implicit none
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:)
  integer :: dummy

  open(unit=20, file=trim(f), status='old', form='unformatted', access='stream')
  read(20) dummy, d1
  close(20)
end subroutine loadfile_1d_d_py

subroutine loadfile_2d_d(f, d1, d2, d3, d4, d5, d6, d7, ns, a)
  implicit none
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:,:)
  double precision, intent(out), optional :: d2(:,:), d3(:,:), d4(:,:), d5(:,:), d6(:,:), d7(:,:)
  logical, intent(in), optional :: ns
  character(*), intent(in), optional :: a
  integer :: frm(2)

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) then
    read(20) frm
    call check_inputdata_size(frm, shape(d1), 'loadfile_2d_d')
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

subroutine loadfile_2d_c(f, d1, d2, d3, d4, d5, d6, d7, ns, a)
  implicit none
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:,:)
  complex(dlc), intent(out), optional :: d2(:,:), d3(:,:), d4(:,:), d5(:,:), d6(:,:), d7(:,:)
  logical, intent(in), optional :: ns
  character(*), intent(in), optional :: a
  integer :: frm(2)

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) then
    read(20) frm
    call check_inputdata_size(frm, shape(d1), 'loadfile_2d_c')
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

subroutine loadfile_2d_d_py(f, d1, d2, d3, d4, d5, d6, d7)
  implicit none
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:,:)
  double precision, intent(out), optional :: d2(:,:), d3(:,:), d4(:,:), d5(:,:), d6(:,:), d7(:,:)
  integer :: dummy

  open(unit=20, file=trim(f), status='old', form='unformatted', access='stream')
  read(20) dummy, d1
  if (present(d2)) read(20) dummy, dummy, d2
  if (present(d3)) read(20) dummy, dummy, d3
  if (present(d4)) read(20) dummy, dummy, d4
  if (present(d5)) read(20) dummy, dummy, d5
  if (present(d6)) read(20) dummy, dummy, d6
  if (present(d7)) read(20) dummy, dummy, d7
  close(20)
end subroutine loadfile_2d_d_py

subroutine loadfile_3d_d(f, d1, d2, d3, d4, d5, d6, d7, ns, a)
  implicit none
  character(*), intent(in) :: f
  double precision, intent(out) :: d1(:,:,:)
  double precision, intent(out), optional :: d2(:,:,:), d3(:,:,:), d4(:,:,:), d5(:,:,:), d6(:,:,:), d7(:,:,:)
  logical, intent(in), optional :: ns
  character(*), intent(in), optional :: a
  integer :: frm(3)

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) then
    read(20) frm
    call check_inputdata_size(frm, shape(d1), 'loadfile_3d_d')
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

subroutine loadfile_3d_c(f, d1, d2, d3, d4, d5, d6, d7, ns, a)
  implicit none
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:,:,:)
  complex(dlc), intent(out), optional :: d2(:,:,:), d3(:,:,:), d4(:,:,:), d5(:,:,:), d6(:,:,:), d7(:,:,:)
  logical, intent(in), optional :: ns
  character(*), intent(in), optional :: a
  integer :: frm(3)

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) then
    read(20) frm
    call check_inputdata_size(frm, shape(d1), 'loadfile_3d_c')
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

subroutine loadfile_4d_c(f, d1, d2, d3, d4, d5, ns, a)
  implicit none
  character(*), intent(in) :: f
  complex(dlc), intent(out) :: d1(:,:,:,:)
  complex(dlc), intent(out), optional :: d2(:,:,:,:), d3(:,:,:,:), d4(:,:,:,:), d5(:,:,:,:)
  logical, intent(in), optional :: ns
  character(*), intent(in), optional :: a
  integer :: frm(4)

  open(unit=20, file=trim(f), status='old', form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) then
    read(20) frm
    call check_inputdata_size(frm, shape(d1), 'loadfile_4d_c')
  end if

  read(20) d1
  if (present(d2)) read(20) d2
  if (present(d3)) read(20) d3
  if (present(d4)) read(20) d4
  if (present(d5)) read(20) d5
  close(20)
end subroutine loadfile_4d_c

subroutine savefile_1d_d(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  double precision, intent(in) :: d1(:)
  double precision, intent(in), optional :: d2(:), d3(:), d4(:), d5(:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) write(20) size(d1)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_1d_d

subroutine savefile_1d_c(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), intent(in) :: d1(:)
  complex(dlc), intent(in), optional :: d2(:), d3(:), d4(:), d5(:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) write(20) size(d1)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_1d_c

subroutine savefile_2d_d(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  double precision, intent(in) :: d1(:,:)
  double precision, intent(in), optional :: d2(:,:), d3(:,:), d4(:,:), d5(:,:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) write(20) size(d1, dim=1), size(d1, dim=2)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_2d_d

subroutine savefile_2d_c(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), intent(in) :: d1(:,:)
  complex(dlc), intent(in), optional :: d2(:,:), d3(:,:), d4(:,:), d5(:,:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) write(20) size(d1, dim=1), size(d1, dim=2)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_2d_c

subroutine savefile_3d_d(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  double precision, intent(in) :: d1(:,:,:)
  double precision, intent(in), optional :: d2(:,:,:), d3(:,:,:), d4(:,:,:), d5(:,:,:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) write(20) size(d1, dim=1), size(d1, dim=2), size(d1, dim=3)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_3d_d

subroutine savefile_3d_c(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), intent(in) :: d1(:,:,:)
  complex(dlc), intent(in), optional :: d2(:,:,:), d3(:,:,:), d4(:,:,:), d5(:,:,:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) write(20) size(d1, dim=1), size(d1, dim=2), size(d1, dim=3)
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_3d_c

subroutine savefile_4d_c(f, d1, d2, d3, d4, d5, ow, ns, a)
  implicit none
  logical, intent(in), optional :: ow, ns
  character(*), intent(in) :: f
  complex(dlc), intent(in) :: d1(:,:,:,:)
  complex(dlc), intent(in), optional :: d2(:,:,:,:), d3(:,:,:,:), d4(:,:,:,:), d5(:,:,:,:)
  character(*), intent(in), optional :: a

  open(unit=20, file=trim(f), status=file_status(ow), form='unformatted', access=access_mode(a))
  if (write_size_header(ns)) then
    write(20) size(d1, dim=1), size(d1, dim=2), size(d1, dim=3), size(d1, dim=4)
  end if
  write(20) d1
  if (present(d2)) write(20) d2
  if (present(d3)) write(20) d3
  if (present(d4)) write(20) d4
  if (present(d5)) write(20) d5
  close(20)
end subroutine savefile_4d_c

end module io
