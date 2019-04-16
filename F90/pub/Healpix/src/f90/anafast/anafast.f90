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
module anamod

  USE healpix_types
  USE alm_tools,  only: map2alm_iterative, alm2cl
  USE fitstools,  only: getsize_fits, input_map, read_par, read_bintab, write_asctab, dump_alms
  USE head_fits,  only: add_card, write_minimal_header
  use misc_utils, only: assert_alloc, fatal_error, wall_clock_time, brag_openmp, string
  USE pix_tools,  only: convert_nest2ring, convert_ring2nest, nside2npix, npix2nside, remove_dipole, vec2ang
  USE paramfile_io, only: paramfile_handle, parse_init, parse_int, &
         parse_string, parse_double, parse_summarize, parse_check_unused, &
         parse_finish, concatnl, scan_directories, get_healpix_data_dir, get_healpix_test_dir, &
         get_healpix_ring_weight_file
  use udgrade_nr, only: udgrade_nest, udgrade_ring
  implicit none

  character(len=*), parameter :: CODE = "ANAFAST"
  character(len=FILENAMELEN)  :: lcode

  private
  public :: ana_sub_s, ana_sub_d, CODE, lcode, check_input_map

  contains

    subroutine ana_sub_s(parafile)
      ! single precision implementation
      integer(I4B),               parameter  :: KMAP  = SP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'ana_sub_inc.f90'
    end subroutine ana_sub_s
    !------------------------------------


    subroutine ana_sub_d(parafile)
      ! double precision implementation
      integer(I4B),               parameter  :: KMAP  = DP  ! precision of map arrays
      character(len=FILENAMELEN), intent(in) :: parafile
      include 'ana_sub_inc.f90'
    end subroutine ana_sub_d
    !------------------------------------


    subroutine check_input_map(mapfile, polarisation)
      ! check out that file contains a valid HEALPIX map
      ! on input polarisation is set to T if one wants to get polarisation data from map,
      ! and on output it will be set to F if that is not possible
      !
      character(len=*)                      :: mapfile
      logical(LGT), intent(inout)           :: polarisation
      !
      integer(I4B) :: nmaps, order_map, nsmax, mlpol, type, polar_fits, polcconv
      integer(I4B) :: npixtot
      character(len=1) :: coordsys

      npixtot = getsize_fits(mapfile, nmaps = nmaps, ordering=order_map, nside=nsmax,&    
           &              mlpol=mlpol, type = type, polarisation = polar_fits, &
           &             coordsys=coordsys, polcconv=polcconv)
      if (nsmax<=0) then
         print*,"Keyword NSIDE not found in FITS header!"
         call fatal_error(code)
      endif
      if (type == 3) npixtot = nside2npix(nsmax) ! cut sky input data set
      if (nsmax/=npix2nside(npixtot)) then
         print 9000,"FITS header keyword NSIDE does not correspond"
         print 9000,"to the size of the map!"
         call fatal_error(code)
      endif

      if (polarisation .and. (nmaps >=3) .and. polar_fits == -1) then
         print 9000,"The input fits file MAY NOT contain polarisation data."
         print 9000,"Proceed at your own risk"
      endif

      if (polarisation .and. (nmaps<3 .or. polar_fits ==0)) then
         print 9000,"The file does NOT contain polarisation maps"
         print 9000,"only the temperature field will be analyzed"
         polarisation = .false.
      endif

      if (polarisation .and. (polcconv == 2)) then
         print 9000,"The input map contains polarized data in the IAU coordinate convention (POLCCONV)"
         print 9000,code//" can not proceed with these data"
         print 9000,"See Healpix primer for details."
         call fatal_error(code)
      endif

      if (polarisation .and. (polcconv == 0)) then
         print 9000,"WARNING: the polarisation coordinate convention (POLCCONV) can not be determined"
         print 9000,"         COSMO will be assumed. See Healpix primer for details."
      endif

      !     --- check ordering scheme ---
      if ((order_map/=1).and.(order_map/=2)) then
         print 9000,"The ordering scheme of the map must be RING or NESTED."
         print 9000,"No ordering specification is given in the FITS-header!"
         call fatal_error(code)
      endif
  
9000  format(a)

      return
    end subroutine check_input_map
      

end module anamod
!====================================================================================

program anafast

  use healpix_types
  use extension, only: getArgument, nArguments
  use anamod,    only: ana_sub_s, ana_sub_d, CODE, lcode, check_input_map
  use misc_utils,only: assert, fatal_error, strlowcase

  integer(i4b)                :: n_args, i
  character(len=FILENAMELEN)  :: arg
  logical(lgt)                :: do_double
  character(len=FILENAMELEN)  :: parafile

  ! count arguments, should be 0, 1 or 2
  n_args = nArguments()
  lcode = trim(strlowcase(CODE))
  call assert(n_args <= 2,' Usage: '//trim(lcode)//' [-s|--single|-d|--double] [parameter_file_name]')

  parafile = ''
  do_double = .false.

  ! parse arguments
  do i=1, n_args
     call getArgument(i, arg)
     arg = trim(adjustl(arg))
     if (arg(1:1) == '-') then
        select case (trim(arg))
        case ('-d', '--double')
           do_double = .true.
        case ('-s', '--single')
           do_double = .false.
        case default
           call fatal_error(CODE//': Invalid argument: '//trim(arg))
        end select
     else
        parafile = arg
     endif
  enddo

  ! start calculations
  if (do_double) then
     call ana_sub_d(parafile)
  else            
     call ana_sub_s(parafile)
  endif

  stop
end program anafast
