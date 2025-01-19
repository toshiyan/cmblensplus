! -*- F90 -*-
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

! ---------------------------------------------------------------------
! Module [gifmod] provides five subroutines: [setcol] to set color
!   tables, [addbar] to add a color bar, [addttl] to add a title
!   string, [imgscl] to scale an image array to the color table,
!   and [gifmap] to write an array as a GIF file.
!
!
! Matthias Bartelmann and Anthony J. Banday, Dec. 1998
! Dec 2010, EH.
!     redefined AS parameter
!     replaced WHERE structure by loops in imgscl
!     used loop to flip A array in gifmap
! 2013-05-07: G95-compatible
! 2019-02-19: simplified C binding (M. Reinecke)
! ---------------------------------------------------------------------

module gifmod
   use iso_c_binding
   use healpix_types
   use misc_utils
   implicit none

!   real(sp), parameter :: as = -1.6375e30, tiny = 1.0e-20
   real(sp),     parameter         :: as = hpx_sbadval, tiny = 1.0e-20
   real(sp),     save                      :: amin,amax
   integer(i4b), save                      :: nc = -1 ! for G95
   integer(i4b), allocatable, dimension(:) :: r,g,b

   ! ------------------------------------------------------------------
   ! pre-defined color tables
   ! ------------------------------------------------------------------

   integer(i4b), dimension(0:5,100) :: r0,g0,b0
   data r0(0,:) /&
        &   0,  2,  5,  7, 10, 12, 15, 18, 20, 23, 25, 28, 30, 33, 36,&
        &  38, 41, 43, 46, 48, 51, 54, 56, 59, 61, 64, 66, 69, 72, 74,&
        &  77, 79, 82, 85, 87, 90, 92, 95, 97,100,103,105,108,110,113,&
        & 115,118,121,123,126,128,131,133,136,139,141,144,146,149,151,&
        & 154,157,159,162,164,167,170,172,175,177,180,182,185,188,190,&
        & 193,195,198,200,203,206,208,211,213,216,218,221,224,226,229,&
        & 231,234,236,239,242,244,247,249,252,255 /
   data g0(0,:) /&
        &   0,  2,  5,  7, 10, 12, 15, 18, 20, 23, 25, 28, 30, 33, 36,&
        &  38, 41, 43, 46, 48, 51, 54, 56, 59, 61, 64, 66, 69, 72, 74,&
        &  77, 79, 82, 85, 87, 90, 92, 95, 97,100,103,105,108,110,113,&
        & 115,118,121,123,126,128,131,133,136,139,141,144,146,149,151,&
        & 154,157,159,162,164,167,170,172,175,177,180,182,185,188,190,&
        & 193,195,198,200,203,206,208,211,213,216,218,221,224,226,229,&
        & 231,234,236,239,242,244,247,249,252,255 /
   data b0(0,:) /&
        &   0,  2,  5,  7, 10, 12, 15, 18, 20, 23, 25, 28, 30, 33, 36,&
        &  38, 41, 43, 46, 48, 51, 54, 56, 59, 61, 64, 66, 69, 72, 74,&
        &  77, 79, 82, 85, 87, 90, 92, 95, 97,100,103,105,108,110,113,&
        & 115,118,121,123,126,128,131,133,136,139,141,144,146,149,151,&
        & 154,157,159,162,164,167,170,172,175,177,180,182,185,188,190,&
        & 193,195,198,200,203,206,208,211,213,216,218,221,224,226,229,&
        & 231,234,236,239,242,244,247,249,252,255 /
   data r0(1,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   5, 13, 25, 33, 45, 57, 65, 77, 85, 98,106,118,130,138,150,&
        & 158,170,179,191,203,211,223,231,243,255 /
   data g0(1,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  2,  7, 12, 15, 20, 23, 28,&
        &  31, 36, 41, 44, 49, 52, 57, 60, 65, 69, 73, 77, 81, 85, 89,&
        &  94, 98,102,106,110,114,119,122,127,130,135,138,143,148,151,&
        & 156,159,164,167,172,177,180,185,188,193,196,201,206,209,214,&
        & 217,222,225,230,235,238,243,246,251,255 /
   data b0(1,:) /&
        &   0,  2,  6,  9, 13, 16, 20, 24, 27, 31, 33, 37, 40, 44, 48,&
        &  51, 55, 58, 62, 65, 69, 73, 75, 80, 82, 86, 89, 93, 97,100,&
        & 104,107,111,115,118,122,124,128,131,135,139,142,146,149,153,&
        & 155,160,164,166,170,173,177,180,184,188,191,195,198,202,204,&
        & 208,212,215,219,222,226,230,233,237,240,244,246,250,255,255,&
        & 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,255 /
   data r0(2,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 12, 30, 48,&
        &  60, 78, 90,108,120,138,156,168,186,198,216,228,243,252,252,&
        & 251,249,248,247,245,242,240,238,236,236,233,231,228,228,227,&
        & 225,222,219,217,216,216,213,211,208,205,204,204,202,199,197,&
        & 196,195,193,190,188,185,184,184,181,179,176,176,175,172,170,&
        & 167,165,164,164,161,158,156,153,152,152,150,151,160,166,176,&
        & 184,194,200,209,218,224,236,243,252,255 /
   data g0(2,:) /&
        &   0, 72, 82, 90,100,108,135,162,180,207,225,252,246,234,216,&
        & 204,186,174,156,144,126,108, 96, 78, 66, 48, 36, 18,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  8, 32, 48, 72,&
        &  88,112,128,149,172,188,212,228,252,255 /
   data b0(2,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   3,  7, 12, 18, 22, 29, 34, 39, 43, 50, 56, 61, 67, 70, 76,&
        &  80, 87, 94, 98,104,108,114,118,125,131,135,141,145,152,156,&
        & 162,167,171,178,183,189,195,199,205,210,216,219,225,232,236,&
        & 243,247,252,255,255,255,255,255,255,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,255 /
   data r0(3,:) /&
        &   0,  2,  7, 10, 14, 17, 21, 26, 28, 33, 36, 40, 43, 47, 52,&
        &  55, 59, 62, 66, 69, 73, 78, 81, 85, 88, 92, 95, 99,104,107,&
        & 111,114,118,123,126,130,133,137,140,144,149,152,156,159,163,&
        & 166,170,175,178,182,185,189,192,197,201,204,208,211,215,218,&
        & 223,227,230,234,237,241,246,249,253,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,255 /
   data g0(3,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  1,  5, 11, 15, 20, 24, 30, 35, 39, 45, 49, 54, 58,&
        &  64, 69, 73, 79, 83, 88, 94, 98,103,107,113,117,122,128,132,&
        & 137,141,147,151,156,162,166,171,175,181,185,190,196,200,205,&
        & 209,215,219,224,230,234,239,243,249,255 /
   data b0(3,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &  11, 19, 31, 39, 51, 62, 70, 82, 90,102,109,121,133,141,153,&
        & 160,172,180,192,204,211,223,231,243,255 /
   data r0(4,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  7,&
        &  22, 45, 67, 82,105,120,135,145,160,175,185,200,201,202,203,&
        & 205,206,207,209,210,211,213,214,215,216,218,219,220,222,223,&
        & 224,225,227,228,229,231,232,233,234,236,237,238,240,241,242,&
        & 243,245,246,247,249,250,251,252,254,255 /
   data g0(4,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3, 12,&
        &  18, 28, 34, 43, 50, 59, 68, 75, 84, 90,100,106,115,125,131,&
        & 140,146,150,150,150,150,150,150,149,147,145,144,142,141,137,&
        & 132,125,117,112,105,100, 81, 68, 50, 31, 18,  0,  4, 11, 16,&
        &  23, 29, 34, 41, 46, 53, 59, 64, 71, 76, 83, 87, 94,101,106,&
        & 113,117,124,129,136,142,147,154,159,166,170,177,184,189,196,&
        & 200,207,212,219,226,230,237,242,249,255 /
   data b0(4,:) /&
        &   0,  4, 10, 14, 20, 25, 31, 37, 41, 47, 52, 58, 62, 68, 75,&
        &  79, 85, 89, 95,100,100,100,100,100,100,100,100,100,100,100,&
        & 100,100, 93, 84, 78, 68, 62, 53, 46, 37, 28, 21, 12,  6,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
   data r0(5,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  7, 19, 27, 39, 47, 59,&
        &  67, 79, 91, 99,111,119,131,139,151,163,171,183,191,203,211,&
        & 223,235,243,255,255,255,255,255,255,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,255,255,255,255,255,241,&
        & 233,219,211,197,184,175,162,153,140,131 /
   data g0(5,:) /&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 11,&
        &  19, 31, 39, 51, 59, 71, 83, 91,103,111,123,131,143,155,163,&
        & 175,183,195,207,215,227,235,247,255,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,&
        & 255,255,255,255,247,235,223,215,203,195,183,175,163,151,143,&
        & 131,123,111,103, 91, 79, 71, 59, 51, 39, 31, 19,  7,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
   data b0(5,:) /&
        & 131,135,147,155,167,175,187,199,207,219,227,239,247,255,255,&
        & 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,&
        & 255,255,255,255,255,255,255,255,255,247,235,227,215,207,195,&
        & 187,175,163,155,143,135,123,115,103, 91, 83, 71, 63, 51, 43,&
        &  31, 19, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,&
        &   0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

   interface setcol
      module procedure setcol_i, setcol_t
   end interface

   private
   public addbar,addttl,setcol,imgscl,gifmap

   interface
     subroutine gifout(x,nx,ny,r,g,b,nc,fn) bind(C, name="gifout")
     use iso_c_binding
     implicit none
     integer (c_int), intent(in) :: x(*), nx, ny, r(*), g(*), b(*), nc
     character(c_char), intent(in) :: fn(*)
     end subroutine
     subroutine gifstr(x,nx,ny,sx,sy,nc,r,g,b,or_,fn) bind(C, name="gifstr")
     use iso_c_binding
     implicit none
     integer (c_int), intent(inout) :: x(*) ! 2020-08-25
     integer (c_int), intent(in) :: nx, ny, sx, sy, nc, r(*), g(*), b(*), or_
     character(c_char), intent(in) :: fn(*)
     end subroutine
   end interface

contains

   ! ------------------------------------------------------------------
   ! gifmap: F90 wrapper around C function gifout;
   !         takes array [a] and produces integer array [x] from it
   !         scaled from 0 to the number of colors minus 1;
   !         output is written to GIF file [fn];
   ! ------------------------------------------------------------------

   subroutine gifmap(a,fn)
      use healpix_types
      implicit none
      integer(i4b) :: nx,ny,j
      integer(i4b), dimension(:,:), intent(in) :: a
      character(len=*), intent(in) :: fn
      integer(i4b), dimension(:,:), allocatable :: arvs

      if (nc <= 0) then
         print *, "error in gifmap: color table undefined"
         call fatal_error
      end if

      nx = size(a,1)
      ny = size(a,2)
      allocate(arvs(nx,ny))
      do j=1, ny
         arvs(:,ny-j+1) = a(:,j)
      enddo
      call gifout(arvs, nx, ny, r, g, b, nc, trim(fn)//C_NULL_CHAR)
      deallocate(arvs)

   end subroutine gifmap

   ! ------------------------------------------------------------------
   ! imgscl: scale real input array [a] to an integer output array [b]
   !         ranging from 3 to the number of colors. The optional mask
   !         [m] specifies what area on the array should be scaled.
   ! ------------------------------------------------------------------

   subroutine imgscl(a,b,m)
      use healpix_types
      implicit none
      real(sp) :: dela
      real(sp), dimension(:,:), intent(in) :: a
      integer(i4b) :: nx,ny
      integer(i4b), dimension(size(a,1),size(a,2)), intent(out) :: b
      logical(lgt), optional, dimension(size(a,1),size(a,2)), intent(in) :: m
      !logical(lgt), allocatable, dimension(:,:) :: n
      real(sp) :: threshold = 1.e-6 * abs(as) ! corrected 2012-08-27
      integer(i4b) :: i,j

      if (present(m)) then
         nx = size(a,1)
         ny = size(a,2)

         amax = -max_sp
         amin =  max_sp
         do j=1, ny
            do i=1, nx
               if (m(i,j) .and. abs(a(i,j)-as) > threshold) then
                  amax = max(a(i,j), amax)
                  amin = min(a(i,j), amin)
               endif
            enddo
         enddo

         dela = max(amax - amin,tiny)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(nx, ny, m, a, threshold, b, amin, amax, dela, nc)&
!$OMP PRIVATE(i, j)
!$OMP DO schedule(dynamic, 1000)
         do j=1, ny
            do i=1, nx
               if (m(i,j)) then
                  if (abs(a(i,j)-as) > threshold) then
                     ! valid pixels
                     b(i,j) = int((a(i,j) - amin) / dela * real(nc - 4)) + 3
                  else
                     ! in map, unobserved pixels
                     b(i,j) = 2
                  endif
               else
                  ! out of map pixels
                  b(i,j) = 1
               endif
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL


      else
         amax = maxval(a)
         amin = minval(a)

         dela = max(amax - amin,tiny)

         b = int((a - amin) / dela * real(nc - 4)) + 3
      end if

   end subroutine imgscl

   ! ------------------------------------------------------------------
   ! setcol: defines arrays [r] (red), [g] (green) and [b] blue
   !         according to predefined color table [i] or according to
   !         input color tables [rc], [gc] and [bc]
   ! ------------------------------------------------------------------

   subroutine setcol_i(i)
      use healpix_types
      implicit none
      integer(i4b), intent(in) :: i

      if (i < 0 .or. i > 5) then
         print *, "error in setcol: undefined color table selected"
         call fatal_error
      end if

      nc = size(r0(i,:)) + 3
      if (allocated(r)) deallocate(r)
      if (allocated(g)) deallocate(g)
      if (allocated(b)) deallocate(b)
      allocate(r(nc),g(nc),b(nc))

      ! reserve three lowest-indexed colors to black, white, and grey

      r(1:3) = (/   0, 255, 175 /)
      g(1:3) = r(1:3)
      b(1:3) = r(1:3)

      ! fill the rest with the predefined colors

      r(4:) = r0(i,:)
      g(4:) = g0(i,:)
      b(4:) = b0(i,:)

   end subroutine setcol_i

   subroutine setcol_t(rc,gc,bc)
      use healpix_types
      implicit none
      integer(i4b), dimension(:), intent(in) :: rc,gc,bc

      nc = size(rc)
      if (size(gc) /= nc .or. size(bc) /= nc) then
         print *, "error in setcol: color tables have different lengths"
         call fatal_error
      end if

      if (allocated(r)) deallocate(r)
      if (allocated(g)) deallocate(g)
      if (allocated(b)) deallocate(b)
      allocate(r(nc),g(nc),b(nc))
      r = rc
      g = gc
      b = bc

   end subroutine setcol_t

   ! ------------------------------------------------------------------
   ! addbar: add a color bar to input array [a], returning the output
   !         in pointer array [b];
   ! ------------------------------------------------------------------

   subroutine addbar(a,b)
      use healpix_types
      implicit none
      integer(i4b) :: nx,ny,my,cx,cy,dx,dy,i,sx,sy
      integer(i4b), dimension(:,:), intent(in) :: a
      integer(i4b), pointer, dimension(:,:) :: b
      real(sp), allocatable, dimension(:,:) :: c
      character(len=12) :: smin,smax

      if (nc <= 0) then
         print *, "error in addbar: color table undefined"
         call fatal_error
      end if

      nx = size(a,1)
      ny = size(a,2)
      my = int(1.1 * real(ny))

      allocate(b(nx,my))
      b = 1 ! minval(a)
      b(:,my-ny+1:my) = a

      cx = int(0.7 * real(nx))
      cy = int(0.3 * real(my - ny))

      dx = (nx - cx) / 2
      dy = (my - ny - cy) / 2

      allocate(c(cx,cy))

      do i=1,cx
        c(i,1) = real(i-1)/real(cx-1)
      end do
      do i = 2,cy
         c(:,i) = c(:,1)
      end do
      c = amin + (amax - amin) * c

      call imgscl(c,b(dx:dx+cx-1,dy:dy+cy-1))

      deallocate(c)

!       write (smin,'(f10.6)') amin
!       write (smax,'(f10.6)') amax
!       if (abs(amin) >= 100.0 .or. abs(amin) <= 0.01)&
!            & write (smin,'(es10.3)') amin
!       if (abs(amax) >= 100.0 .or. abs(amax) <= 0.01)&
!            & write (smax,'(es10.3)') amax
      write (smin,'(1pg11.3)') amin
      write (smax,'(1pg11.3)') amax

      if (amin > 0.0) smin = '+' // trim(adjustl(smin))
      if (amax > 0.0) smax = '+' // trim(adjustl(smax))

      sx = dx
      sy = my - dy - cy / 2
      call addstr(b,sx,sy,-1,smin)

      sx = dx + cx - 1
      sy = my - dy - cy / 2
      call addstr(b,sx,sy, 1,smax)

   end subroutine addbar

   ! ------------------------------------------------------------------
   ! addttl: add title string [st] to map [a], return new image [b]
   ! ------------------------------------------------------------------

   subroutine addttl(a,b,st)
      use healpix_types
      implicit none
      integer(i4b) :: nx,ny,my,sx,sy
      integer(i4b), dimension(:,:), intent(in) :: a
      integer(i4b), pointer, dimension(:,:) :: b
      character(len=*), intent(in) :: st

      nx = size(a,1)
      ny = size(a,2)
      my = int(1.1 * real(ny))

      allocate(b(nx,my))
      b = 1
      b(:,1:ny) = a

      sx = nx / 2
      sy = (my - ny) / 2

      call addstr(b,sx,sy, 0,st)

   end subroutine addttl

   ! ------------------------------------------------------------------
   ! addstr: used by addbar to add labels to the ends of color bar 
   !          and by addttl to add title; 
   !         [a] is the (modified) image array, [i] and [j] are the
   !         array positions where the string [s] should end if [or]
   !         is negative, or start if [or] is positive
   ! ------------------------------------------------------------------

   subroutine addstr(a,i,j,or,s)
      use healpix_types
      implicit none
      integer(i4b), intent(inout), dimension(:,:) :: a
      integer(i4b), intent(in) :: i,j
      integer(i4b), intent(in) :: or
      character(len=*), intent(in) :: s
      integer(i4b) :: nx,ny

      nx = size(a,1)
      ny = size(a,2)

      call gifstr(a(:,ny:1:-1),nx,ny,i,j,nc,r,g,b,or,trim(s)//C_NULL_CHAR)

   end subroutine addstr

end module gifmod
