!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.sourceforge.net
!
!-----------------------------------------------------------------------------
! template for routine SP/DP overloading for module mpi_alm_tools

! K M A P   : map kind                 either SP or DP
! K A L M C : alm kind (complex)       either SPC or DPC
! K A L M   : alm related and cl kind  either SP or DP
!
! K L O A D : suffix of routine name, to be replaced by either s (SP) or d (DP)
!
!


  subroutine mpi_map2alm_simple_KLOAD(comm_in, map, alms, zbounds_in, w8ring)
    implicit none
    
    integer(i4b),                        intent(in)           :: comm_in
    real(KMAP),     dimension(0:,1:),    intent(in)           :: map
    complex(KALMC), dimension(1:,0:,0:), intent(out)          :: alms
    real(dp),       dimension(2),        intent(in), optional :: zbounds_in
    real(dp),       dimension(1:,1:),    intent(in), optional :: w8ring

    integer(i4b) :: myid_int, nmaps_int, npix_int, nsmax_in, nlmax_in, nmmax_in
    integer(i4b) :: precompute_plms_in
    logical(lgt) :: polarization_in
    real(dp),              dimension(2)   :: zbounds_int
    real(dp), allocatable, dimension(:,:) :: w8ring_int

    call mpi_comm_rank(comm_in, myid_int, ierr)

    if (myid_int == root) then

       npix_int   = size(map(:,1))
       nmaps_int  = size(map(0,:))
       nsmax_in   = nint(sqrt(real(npix_int,sp)/12.))
       nlmax_in   = size(alms(1,:,0))-1
       nmmax_in   = nlmax_in
       if (nmaps_int == 3) then
          polarization_in = .true.
       else
          polarization_in = .false.
       end if
       precompute_plms_in = 0
       
       allocate(w8ring_int(1:2*nsmax_in,nmaps_int))
       if (present(w8ring)) then
          w8ring_int = w8ring
       else
          w8ring_int = 1.d0
       end if
       
       if (present(zbounds_in)) then
          zbounds_int = zbounds_in
       else
          zbounds_int = 0.d0
       end if

       call mpi_initialize_alm_tools(comm_in, nsmax_in, nlmax_in, &
            & nmmax_in, zbounds_int, polarization_in, precompute_plms_in, &
            & w8ring_int)

       call mpi_map2alm(map, alms)

       deallocate(w8ring_int)

    else

       call mpi_initialize_alm_tools(comm_in)

       call mpi_map2alm_slave

    end if

    call mpi_cleanup_alm_tools

  end subroutine mpi_map2alm_simple_KLOAD


  subroutine mpi_alm2map_simple_KLOAD(comm_in, alms, map)
    implicit none
    
    integer(i4b),                        intent(in)  :: comm_in
    complex(KALMC), dimension(1:,0:,0:), intent(in)  :: alms
    real(KMAP),     dimension(0:,1:),    intent(out) :: map

    integer(i4b) :: myid_int, npix_int, nmaps_int, nsmax_in, nlmax_in, nmmax_in
    integer(i4b) :: precompute_plms_in
    logical(lgt) :: polarization_in
    real(dp),     allocatable, dimension(:,:)  :: w8ring_int
    real(dp),     dimension(2)                 :: zbounds_int

    call mpi_comm_rank(comm_in, myid_int, ierr)

    if (myid_int == root) then

       nmaps_int   = size(map(0,:))
       npix_int    = size(map(:,1))
       nsmax_in    = nint(sqrt(real(npix_int,sp)/12.))
       nlmax_in    = size(alms(1,:,0))-1
       nmmax_in    = nlmax_in
       zbounds_int = 0.d0
       if (nmaps_int == 3) then
          polarization_in = .true.
       else
          polarization_in = .false.
       end if
       precompute_plms_in = 0
       
       allocate(w8ring_int(1:2*nsmax_in,nmaps_int))
       w8ring_int = 1.d0

       call mpi_initialize_alm_tools(comm_in, nsmax_in, nlmax_in, &
            & nmmax_in, zbounds_int, polarization_in, precompute_plms_in, &
            & w8ring_int)

       call mpi_alm2map(alms, map)

       deallocate(w8ring_int)

    else

       call mpi_initialize_alm_tools(comm_in)

       call mpi_alm2map_slave

    end if

    call mpi_cleanup_alm_tools

  end subroutine mpi_alm2map_simple_KLOAD



  ! ***************************************************************
  !            Wrappers for the computational routines
  ! ***************************************************************

  subroutine mpi_map2alm_KLOAD(map, alms)
    implicit none

    real(KMAP),     dimension(0:,1:),    intent(in)  :: map
    complex(KALMC), dimension(1:,0:,0:), intent(out) :: alms

    call distribute_map(map)

    local_alms = cmplx(0.d0,0.d0)
    alms_work  = cmplx(0.d0,0.d0)
    
    if (polarization) then
       ! Polarization routines
       if (precompute_plms == 1) then
          call mpi_map2alm_pol_pre1
       else if (precompute_plms == 2) then
          call mpi_map2alm_pol_pre2
       else
          call mpi_map2alm_pol
       end if
    else
       ! Temperature routines
       if (precompute_plms == 1 .or. precompute_plms == 2) then
          call mpi_map2alm_sc_pre
       else
          call mpi_map2alm_sc
       end if
    end if

    call mpi_reduce(local_alms, alms_work, size(local_alms), &
         & MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm, ierr) 

    alms = alms_work

  end subroutine mpi_map2alm_KLOAD


  subroutine mpi_alm2map_KLOAD(alms, map)
    implicit none

    complex(KALMC), dimension(1:,0:,0:), intent(in)  :: alms
    real(KMAP),     dimension(0:,1:),    intent(out) :: map

    local_alms = alms

    call mpi_bcast(local_alms, size(local_alms), MPI_DOUBLE_COMPLEX, root, comm, ierr)

    if (polarization) then
       ! Polarization routines
       if (precompute_plms == 1) then
          call mpi_alm2map_pol_pre1
       else if (precompute_plms == 2) then
          call mpi_alm2map_pol_pre2
       else
          call mpi_alm2map_pol
       end if
    else
       ! Temperature routines
       if (precompute_plms == 1 .or. precompute_plms == 2) then
          call mpi_alm2map_sc_pre
       else
          call mpi_alm2map_sc
       end if
    end if

    call gather_map(map)

  end subroutine mpi_alm2map_KLOAD



  ! ***************************************************************
  !                Internal distribution routines
  ! ***************************************************************


  subroutine distribute_map_KLOAD(map)
    ! Distributes the map data
    implicit none

    real(KMAP), dimension(0:,1:), intent(in)  :: map

    integer(i4b) :: i, seg_size, ierr
    integer(i4b), dimension(MPI_STATUS_SIZE) :: status
 
    !Process 0 distributes the map data to the other proc's
    do i = 1, numprocs-1
       
       if (num_segments(i) == 2) then
          
          seg_size = segments(0,1,i)-segments(0,0,i)+1

          buffer(0:seg_size-1,1:nmaps) = &
               & real(map(segments(0,0,i):segments(0,1,i),1:nmaps),dp)
          buffer(seg_size:2*seg_size-1,1:nmaps) = &
               & real(map(segments(1,0,i):segments(1,1,i),1:nmaps),dp)
          call mpi_send(buffer(0:2*seg_size-1,1:nmaps), &
               & 2*seg_size*nmaps, MPI_DOUBLE_PRECISION, i, 90, comm, ierr)
          
       else !Equatorial section -> only 1 segment
          
          seg_size = segments(0,1,i)-segments(0,0,i)+1

          buffer(0:seg_size-1,1:nmaps) = real(map(segments(0,0,i):segments(0,1,i),1:nmaps),dp)
          call mpi_send(buffer(0:seg_size-1,1:nmaps), &
               & seg_size*nmaps, MPI_DOUBLE_PRECISION, i, 90, comm, ierr)
          
       end if
       
    end do
    
    if (num_segments(0) == 2) then

       seg_size = segments(0,1,0)-segments(0,0,0)+1

       local_map(0:seg_size-1,1:nmaps) = &
            & real(map(segments(0,0,0):segments(0,1,0),1:nmaps),dp)
       local_map(seg_size:2*seg_size-1,1:nmaps) = &
            & real(map(segments(1,0,0):segments(1,1,0),1:nmaps),dp)

    else

       seg_size = segments(0,1,0)-segments(0,0,0)+1
       local_map(0:seg_size-1,1:nmaps) = &
            & real(map(segments(0,0,0):segments(0,1,0),1:nmaps),dp)

    end if

  end subroutine distribute_map_KLOAD


  subroutine gather_map_KLOAD(map)
    ! Gathers the map data
    implicit none

    real(KMAP), dimension(0:,1:), intent(out) :: map

    integer(i4b) :: i, seg_size, ierr, maxsize
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: status

    !Process 0 recieves data from the other proc's
    do i=1,numprocs-1
          
       if (num_segments(i) == 2) then
          
          seg_size = segments(0,1,i)-segments(0,0,i)+1

          call mpi_recv(buffer(0:2*seg_size-1,1:nmaps), &
               & 2*seg_size*nmaps, MPI_DOUBLE_PRECISION, i, 91, comm, status, ierr)
          map(segments(0,0,i):segments(0,1,i),1:nmaps) = &
               & real(buffer(0:seg_size-1,1:nmaps),KMAP)
          map(segments(1,0,i):segments(1,1,i),1:nmaps) = &
               & real(buffer(seg_size:2*seg_size-1,1:nmaps),KMAP)
          
       else !Equatorial section -> only 1 segment
          
          seg_size = segments(0,1,i)-segments(0,0,i)+1

          call mpi_recv(buffer(0:seg_size-1,1:nmaps), &
               & seg_size*nmaps, MPI_DOUBLE_PRECISION, i, 91, comm, status, ierr)
          map(segments(0,0,i):segments(0,1,i),1:nmaps) = real(buffer(0:seg_size-1,1:nmaps),KMAP)
          
       end if
       
    end do
    
    if (num_segments(0) == 2) then
       
       seg_size = segments(0,1,0)-segments(0,0,0)+1

       map(segments(0,0,0):segments(0,1,0),1:nmaps) = &
            & real(local_map(0:seg_size-1,1:nmaps),KMAP)       
       map(segments(1,0,0):segments(1,1,0),1:nmaps) = &         
            & real(local_map(seg_size:2*seg_size-1,1:nmaps),KMAP)
       
    else
       
       seg_size = segment(0,1)-segment(0,0)+1
       map(segment(0,0):segment(0,1),1:nmaps) = real(local_map(0:seg_size-1,1:nmaps),KMAP)
       
    end if

  end subroutine gather_map_KLOAD


