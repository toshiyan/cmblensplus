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

MODULE healpix_types
  ! This module sets the types used in the Fortran 90 modules
  ! of the HEALPIX distribution and follows the example of Numerical Recipes
  !
  ! Benjamin D. Wandelt October 1997
  ! Eric Hivon June 1998
  ! Eric Hivon Oct  2001, edited to be compatible with 'F' compiler
  ! Eric Hivon July 2002, addition of i8b, i2b, i1b
  !                       addition of max_i8b, max_i2b and max_i1b
  !            Jan 2005, explicit form of max_i1b because of ifc 8.1.021
  !            June 2005, redefine i8b as 16 digit integer because of Nec f90 compiler
  !            Mars 2008: i8b same as i4b on machines not supporting 64 bits (NO64BITS flag set)
  !            Feb  2009: introduce healpix_version
  !
  character(len=*), PARAMETER, public :: healpix_version = '3.31'
  INTEGER, PARAMETER, public :: i4b = SELECTED_INT_KIND(9)
#ifdef NO64BITS
  INTEGER, PARAMETER, public :: i8b = i4b
#else
  INTEGER, PARAMETER, public :: i8b = SELECTED_INT_KIND(16)
#endif
  INTEGER, PARAMETER, public :: i2b = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER, public :: i1b = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER, public :: sp  = SELECTED_REAL_KIND(5,30)
  INTEGER, PARAMETER, public :: dp  = SELECTED_REAL_KIND(12,200)
  INTEGER, PARAMETER, public :: lgt = KIND(.TRUE.)
  INTEGER, PARAMETER, public :: spc = KIND((1.0_sp, 1.0_sp))
  INTEGER, PARAMETER, public :: dpc = KIND((1.0_dp, 1.0_dp))
  !
  INTEGER(I8B),  PARAMETER, public :: max_i8b = HUGE(1_i8b)
  INTEGER,       PARAMETER, public :: max_i4b = HUGE(1_i4b)
  INTEGER,       PARAMETER, public :: max_i2b = HUGE(1_i2b)
  INTEGER,       PARAMETER, public :: max_i1b = 127
  REAL(kind=sp), PARAMETER, public :: max_sp  = HUGE(1.0_sp)
  REAL(kind=dp), PARAMETER, public :: max_dp  = HUGE(1.0_dp)

  ! Numerical Constant (Double precision)
  REAL(kind=dp), PARAMETER, public :: QUARTPI=0.785398163397448309615660845819875721049_dp
  REAL(kind=dp), PARAMETER, public :: HALFPI= 1.570796326794896619231321691639751442099_dp
  REAL(kind=dp), PARAMETER, public :: PI    = 3.141592653589793238462643383279502884197_dp
  REAL(kind=dp), PARAMETER, public :: TWOPI = 6.283185307179586476925286766559005768394_dp
  REAL(kind=dp), PARAMETER, public :: FOURPI=12.56637061435917295385057353311801153679_dp
  REAL(kind=dp), PARAMETER, public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
  REAL(kind=dp), PARAMETER, public :: EULER = 0.5772156649015328606065120900824024310422_dp
  REAL(kind=dp), PARAMETER, public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
  REAL(kind=dp), PARAMETER, public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_dp

  real(kind=DP), parameter, public :: RAD2DEG = 180.0_DP / PI
  real(kind=DP), parameter, public :: DEG2RAD = PI / 180.0_DP
  real(kind=SP), parameter, public :: hpx_sbadval = -1.6375e30_sp
  real(kind=DP), parameter, public :: hpx_dbadval = -1.6375e30_dp

  ! Maximum length of filenames
  integer, parameter :: filenamelen = 1024

! ! ---- Normalisation and convention ----
  real(kind=dp), parameter, public ::  KvS = 1.0_dp ! 1.0 : CMBFAST (Healpix 1.2)

END MODULE healpix_types


module healpix_fft
  use healpix_types
  use extension, only: exit_with_status
  implicit none
  private

  ! module for FFT operations
  ! edited for compatibility with Absoft Pro, G95 and GFORTRAN
  ! edited Sept 6, 2005 for GFORTRAN
  !
  public real_fft2, complex_fft2, make_fft2_plan, destroy_fft2_plan, planck_fft2_plan, init_fft2_plan, complex_fft, real_fft

  ! routines of Healpix 2.0, compatible with FFTW 2.x

  logical, parameter, public :: fft2_forward=.false., fft2_backward=.true.

  type planck_fft2_plan
     logical      :: direction
     integer(i4b) :: length
  end type planck_fft2_plan

  interface complex_fft2
     module procedure d_c_complex_fft2, d_r_complex_fft2
  end interface

  interface real_fft2
     module procedure s_real_fft2, d_real_fft2
  end interface

  ! routines of Healpix 1.2

  interface complex_fft
     module procedure complex_fft_orig, complex_fft_alt
  end interface

  interface real_fft
     module procedure s_real_fft, d_real_fft
  end interface

contains

  subroutine fft_gpd (data,nn,backward,onlyreal)
    real(dp), intent(inout) :: data(*)
!    integer(i4b), intent(in) :: nn(:)
    integer, intent(in) :: nn(:)
    logical, intent(in) :: backward, onlyreal

!  real(dp) :: work(2*maxval(nn)) ! not valid with Absoft
   real(dp) :: work(2*nn(1)) ! assumes 1 dimension
  real(dp) :: difr, difi, oldsr, oldsi, sumr, sumi, sinth
  real(dp) :: twowr, t2r,t2i, t3r, t3i, t4r, t4i
  real(dp) :: theta
  real(dp) :: tempr, tempi
  real(dp) :: u1r, u1i, u2r, u2i, u3r, u3i, u4r, u4i
  real(dp) :: wr, wi, w2r, w2i, w3r, w3i, wstpr, wstpi
  integer(i4b) :: i, imin, i1, i2, i3, imax, i1max, i2max, idiv, idim

  integer(i4b) :: ifact(32), iff, ifp1, ifp2, ipar, i1rng, icase, irem, iquot
  integer(i4b) :: j, jmin, j1, j2, j3, j1min, j2min, j2stp
  integer(i4b) :: jmax, j1max, j2max, j3max, j1cnj, j1rng, j1rg2
  integer(i4b) :: k, kconj, krang, kmin, k1, k2, k3, k4, kdif,kstep
  integer(i4b) :: l, lmax
  integer(i4b) :: m, mmax
  integer(i4b) :: n, nhalf, np0, np1, np2, np1hf, np2hf, non2, non2t, ntwo
  integer(i4b) :: ntot, nprev
  integer(i4b) :: ndim

!=======================================================================
  ndim = size(nn)

  nprev = 0
  np0 = 0
  wr = 0d0
  wi = 0d0
  w2r = 0d0
  w2i = 0d0
  w3r = 0d0
  w3i = 0d0
  wstpr = 0d0
  wstpi = 0d0

  if (ndim<1) return
  ntot=2
  do idim=1,ndim
    if (nn(idim)<=0) return
    if (2*nn(idim) > size(work)) call exit_with_status(1,"FFT_GDP: work array too small")
    ntot=ntot*nn(idim)
  end do

!---
!--- MAIN LOOP FOR EACH DIMENSION
!---
  np1 = 2
  dimloop: do idim=1, ndim
    n = nn(idim)
    np2 = np1*n
    if (n==1) goto 900
!---
!--- FACTOR N
!---
    m = n
    ntwo = np1
    iff = 1
    idiv = 2
    do
      iquot = m/idiv
      irem = m-idiv*iquot
      if (iquot<idiv) then
        if (irem/=0) then
          ifact(iff) = m
        else
          ntwo = ntwo+ntwo
        endif
        goto 70
      endif
      if (irem/=0) exit
      ntwo = ntwo+ntwo
      m = iquot
    end do

    idiv = 3
    do
      iquot = m/idiv
      irem = m-idiv*iquot
      if (iquot<idiv) then
        ifact(iff) = m
        goto 70
      endif
      if (irem==0) then
        ifact (iff) = idiv
        iff = iff+1
        m = iquot
      else
        idiv = idiv+2
      endif
    end do

!---
!--- SEPARATE FOUR CASES--
!1. COMPLEX TRANSFORM OR REAL TRANSFORM FOR THE 4TH, 5TH,ETC.
!   DIMENSIONS.
!2. REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION.  METHOD--
!   TRANSFORM HALF THE DATA, SUPPLYING THE OTHER HALF BY CON-
!   JUGATE SYMMETRY.
!3. REAL TRANSFORM FOR THE 1ST DIMENSION, N ODD.  METHOD--
!   TRANSFORM HALF THE DATA AT EACH STAGE, SUPPLYING THE OTHER
!   HALF BY CONJUGATE SYMMETRY.
!4. REAL TRANSFORM FOR THE 1ST DIMENSION, N EVEN.  METHOD--
!   TRANSFORM A COMPLEX ARRAY OF LENGTH N/2 WHOSE REAL PARTS
!   ARE THE EVEN NUMBERED REAL VALUES AND WHOSE IMAGINARY PARTS
!   ARE THE ODD NUMBERED REAL VALUES.  SEPARATE AND SUPPLY
!   THE SECOND HALF BY CONJUGATE SYMMETRY.
!---
70  non2 = np1*(np2/ntwo)
    if ((idim>=4) .or. (.not. onlyreal)) then
      icase = 1
    elseif (idim>1) then
      icase = 2
    elseif (ntwo<=np1) then
      icase = 3
    else
      icase = 4
      ntwo = ntwo/2
      n = n/2
      np2 = np2/2
      ntot = ntot/2
      i = 3
      do j=2,ntot
        data(j) = data(i)
        i = i+2
      end do
    endif
    i1rng = np1
    if (icase==2) i1rng = np0*(1+nprev/2)
!---
!--- SHUFFLE ON THE FACTORS OF TWO IN N.  AS THE SHUFFLING
!--- CAN BE DONE BY SIMPLE INTERCHANGE, NO WORKING ARRAY IS NEEDED
!---
    if (ntwo<=np1) goto 600
    np2hf = np2/2
    j = 1
    do i2=1,np2,non2
      if (j<i2) then
        i1max = i2+non2-2
        do i1=i2,i1max,2
          do i3=i1,ntot,np2
            j3 = j+i3-i2
            tempr = data(i3)
            tempi = data(i3+1)
            data(i3) = data(j3)
            data(i3+1) = data(j3+1)
            data(j3) = tempr
            data(j3+1) = tempi
          end do
        end do
      endif
      m = np2hf
      do
        if (j<=m) exit
        j = j-m
        m = m/2
        if (m<non2) exit
      end do
      j = j+m
    end do
!---
!--- MAIN LOOP FOR FACTORS OF TWO.  PERFORM FOURIER TRANSFORMS OF
!--- LENGTH FOUR, WITH ONE OF LENGTH TWO IF NEEDED.  THE TWIDDLE FACTOR
!--- W=EXP(ISIGN*2*PI*SQRT(-1)*M/(4*MMAX)).  CHECK FOR W=ISIGN*SQRT(-1)
!--- AND REPEAT FOR W=ISIGN*SQRT(-1)*CONJUGATE(W).
!---
    non2t = non2+non2
    ipar = ntwo/np1
    do while (ipar>2)
      ipar = ipar/4
    end do
    if (ipar==2) then
      do i1=1,i1rng,2
        do j3=i1,non2,np1
          do k1=j3,ntot,non2t
            k2 = k1+non2
            tempr = data(k2)
            tempi = data(k2+1)
            data(k2) = data(k1)-tempr
            data(k2+1) = data(k1+1)-tempi
            data(k1) = data(k1)+tempr
            data(k1+1) = data(k1+1)+tempi
          end do
        end do
      end do
    endif
    mmax = non2

    do while (mmax<np2hf)
      lmax = max(non2t,mmax/2)
      if (mmax>non2) then
        theta = -twopi*real(non2,kind=dp)/real(4*mmax,kind=dp)
        if (backward) theta = -theta
        wr = cos(theta)
        wi = sin(theta)
        wstpr = -2d0*wi*wi
        wstpi = 2d0*wr*wi
      endif
      do l=non2,lmax,non2t
        m = l
        if (mmax<=non2) goto 420
410     w2r = wr*wr-wi*wi
        w2i = 2d0*wr*wi
        w3r = w2r*wr-w2i*wi
        w3i = w2r*wi+w2i*wr
420     do i1=1,i1rng,2
          do j3=i1,non2,np1
            kmin = j3+ipar*m
            if (mmax<=non2) kmin = j3
            kdif = ipar*mmax
            do
              kstep = 4*kdif
              do k1=kmin,ntot,kstep
                k2 = k1+kdif
                k3 = k2+kdif
                k4 = k3+kdif
                if (mmax<=non2) then
                  u1r = data(k1)+data(k2)
                  u1i = data(k1+1)+data(k2+1)
                  u2r = data(k3)+data(k4)
                  u2i = data(k3+1)+data(k4+1)
                  u3r = data(k1)-data(k2)
                  u3i = data(k1+1)-data(k2+1)
                  if (.not. backward) then
                    u4r = data(k3+1)-data(k4+1)
                    u4i = data(k4)-data(k3)
                  else
                    u4r = data(k4+1)-data(k3+1)
                    u4i = data(k3)-data(k4)
                  endif
                else
                  t2r = w2r*data(k2)-w2i*data(k2+1)
                  t2i = w2r*data(k2+1)+w2i*data(k2)
                  t3r = wr*data(k3)-wi*data(k3+1)
                  t3i = wr*data(k3+1)+wi*data(k3)
                  t4r = w3r*data(k4)-w3i*data(k4+1)
                  t4i = w3r*data(k4+1)+w3i*data(k4)
                  u1r = data(k1)+t2r
                  u1i = data(k1+1)+t2i
                  u2r = t3r+t4r
                  u2i = t3i+t4i
                  u3r = data(k1)-t2r
                  u3i = data(k1+1)-t2i
                  if (.not. backward) then
                    u4r = t3i-t4i
                    u4i = t4r-t3r
                  else
                    u4r = t4i-t3i
                    u4i = t3r-t4r
                  endif
                endif
                data(k1) = u1r+u2r
                data(k1+1) = u1i+u2i
                data(k2) = u3r+u4r
                data(k2+1) = u3i+u4i
                data(k3) = u1r-u2r
                data(k3+1) = u1i-u2i
                data(k4) = u3r-u4r
                data(k4+1) = u3i-u4i
              end do
              kmin = 4*(kmin-j3)+j3
              kdif = kstep
              if (kdif>=np2) exit
            end do
          end do
        end do
        m = mmax-m
        if (.not. backward) then
          tempr = wr
          wr = -wi
          wi = -tempr
        else
          tempr = wr
          wr = wi
          wi = tempr
        endif
        if (m>lmax) goto 410
        tempr = wr
        wr = wr*wstpr-wi*wstpi+wr
        wi = wi*wstpr+tempr*wstpi+wi
      end do
      ipar = 3-ipar
      mmax = mmax+mmax
    end do
!---
!--- MAIN LOOP FOR FACTORS NOT EQUAL TO TWO.  APPLY THE TWIDDLE FACTOR
!--- W=EXP(ISIGN*2*PI*SQRT(-1)*(J2-1)*(J1-J2)/(NP2*IFP1)), THEN
!--- PERFORM A FOURIER TRANSFORM OF LENGTH IFACT(IFf), MAKING USE OF
!--- CONJUGATE SYMMETRIES.
!---
600 if (ntwo>=np2) goto 700
    ifp1 = non2
    iff = 1
    np1hf = np1/2
610 ifp2 = ifp1/ifact(iff)
    j1rng = np2
    if (icase==3) then
      j1rng = (np2+ifp1)/2
      j2stp = np2/ifact(iff)
      j1rg2 = (j2stp+ifp2)/2
    endif
    j2min = 1+ifp2
    if (ifp1<np2) then
      do j2=j2min,ifp1,ifp2
        theta = -twopi*real(j2-1,kind=dp)/real(np2,kind=dp)
        if (backward) theta = -theta
        sinth = sin(theta/2d0)
        wstpr = -2d0*sinth*sinth
        wstpi = sin(theta)
        wr = wstpr+1d0
        wi = wstpi
        j1min = j2+ifp1
        do j1=j1min,j1rng,ifp1
          i1max = j1+i1rng-2
          do i1=j1,i1max,2
            do i3=i1,ntot,np2
              j3max = i3+ifp2-np1
              do j3=i3,j3max,np1
                tempr = data(j3)
                data(j3) = data(j3)*wr-data(j3+1)*wi
                data(j3+1) = tempr*wi+data(j3+1)*wr
              end do
            end do
          end do
          tempr = wr
          wr = wr*wstpr-wi*wstpi+wr
          wi = tempr*wstpi+wi*wstpr+wi
        end do
      end do
    endif
    theta = -twopi/real(ifact(iff),kind=dp)
    if (backward) theta = -theta
    sinth = sin(theta/2d0)
    wstpr = -2d0*sinth*sinth
    wstpi = sin(theta)
    kstep = 2*n/ifact(iff)
    krang = kstep*(ifact(iff)/2)+1
    do i1=1,i1rng,2
      do i3=i1,ntot,np2
        do kmin=1,krang,kstep
          j1max = i3+j1rng-ifp1
          do j1=i3,j1max,ifp1
            j3max = j1+ifp2-np1
            do j3=j1,j3max,np1
              j2max = j3+ifp1-ifp2
              k = kmin+(j3-j1+(j1-i3)/ifact(iff))/np1hf
              if (kmin<=1) then
                sumr = 0d0
                sumi = 0d0
                do j2=j3,j2max,ifp2
                  sumr = sumr+data(j2)
                  sumi = sumi+data(j2+1)
                end do
                work(k) = sumr
                work(k+1) = sumi
              else
                kconj = k+2*(n-kmin+1)
                j2 = j2max
                sumr = data(j2)
                sumi = data(j2+1)
                oldsr = 0d0
                oldsi = 0d0
                j2 = j2-ifp2
                do
                  tempr = sumr
                  tempi = sumi
                  sumr = twowr*sumr-oldsr+data(j2)
                  sumi = twowr*sumi-oldsi+data(j2+1)
                  oldsr = tempr
                  oldsi = tempi
                  j2 = j2-ifp2
                  if (j2<=j3) exit
                end do
                tempr = wr*sumr-oldsr+data(j2)
                tempi = wi*sumi
                work(k) = tempr-tempi
                work(kconj) = tempr+tempi
                tempr = wr*sumi-oldsi+data(j2+1)
                tempi = wi*sumr
                work(k+1) = tempr+tempi
                work(kconj+1) = tempr-tempi
              endif
            end do
          end do
          if (kmin<=1) then
            wr = wstpr+1d0
            wi = wstpi
          else
            tempr = wr
            wr = wr*wstpr-wi*wstpi+wr
            wi = tempr*wstpi+wi*wstpr+wi
          endif
          twowr=wr+wr
        end do
        if ((icase/=3) .or. (ifp1>=np2)) then
          k = 1
          i2max = i3+np2-np1
          do i2=i3,i2max,np1
            data(i2) = work(k)
            data(i2+1) = work(k+1)
            k = k+2
          end do
        else
!---
!--- COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N ODD, BY CON-
!--- JUGATE SYMMETRIES AT EACH STAGE.
!---
          j3max = i3+ifp2-np1
          do j3=i3,j3max,np1
            j2max = j3+np2-j2stp
            do j2=j3,j2max,j2stp
              j1max = j2+j1rg2-ifp2
              j1cnj = j3+j2max+j2stp-j2
              do j1=j2,j1max,ifp2
                k = 1+j1-i3
                data(j1) = work(k)
                data(j1+1) = work(k+1)
                if (j1>j2) then
                  data(j1cnj) = work(k)
                  data(j1cnj+1)= -work(k+1)
                endif
                j1cnj = j1cnj-ifp2
              end do
            end do
          end do
        endif
      end do
    end do
    iff = iff+1
    ifp1 = ifp2
    if (ifp1>np1) goto 610
!---
!--- COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N EVEN, BY CON-
!--- JUGATE SYMMETRIES.
!---

700 if ((icase==1) .or. (icase==3)) goto 900
    if (icase==2) goto 800
    nhalf = n
    n = n+n
    theta = -twopi/real(n,kind=dp)
    if (backward) theta = -theta
    sinth = sin (theta/2d0)
    wstpr = -2d0*sinth*sinth
    wstpi = sin(theta)
    wr = wstpr+1d0
    wi = wstpi
    imin = 3
    jmin = 2*nhalf-1
    do while (imin<jmin)
      j = jmin
      do i=imin,ntot,np2
        sumr = (data(i)+data(j))/2d0
        sumi = (data(i+1)+data(j+1))/2d0
        difr = (data(i)-data(j))/2d0
        difi = (data(i+1)-data(j+1))/2d0
        tempr = wr*sumi+wi*difr
        tempi = wi*sumi-wr*difr
        data(i) = sumr+tempr
        data(i+1) = difi+tempi
        data(j) = sumr-tempr
        data(j+1) = -difi+tempi
        j = j+np2
      end do
      imin = imin+2
      jmin = jmin-2
      tempr = wr
      wr = wr*wstpr-wi*wstpi+wr
      wi = tempr*wstpi+wi*wstpr+wi
    end do
    if ((imin<=jmin) .and. (.not. backward)) then
      do i=imin,ntot,np2
        data(i+1) = -data(i+1)
      end do
    endif
    np2 = np2+np2
    ntot = ntot+ntot
    j = ntot+1
    imax = ntot/2+1
    do
      imin = imax-2*nhalf
      i = imin

      do
        i = i+2
        j = j-2
        if (i>=imax) exit
        data(j) = data(i)
        data(j+1) = -data(i+1)
      end do

      data(j) = data(imin)-data(imin+1)
      data(j+1) = 0d0

      if (i>=j) exit

      do
        i = i-2
        j = j-2
        if (i<=imin) exit
        data(j) = data(i)
        data(j+1) = data(i+1)
      end do
      data(j) = data(imin)+data(imin+1)
      data(j+1) = 0d0
      imax = imin
    end do
    data(1) = data(1)+data(2)
    data(2) = 0d0
    goto 900
!---
!--- COMPLETE A REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION BY
!--- CONJUGATE SYMMETRIES.
!---
800 if (i1rng>=np1) goto 900
    do i3=1,ntot,np2
      i2max = i3+np2-np1
      do i2=i3,i2max,np1
        imin = i2+i1rng
        imax = i2+np1-2
        jmax = 2*i3+np1-imin
        if (i2>i3) jmax = jmax+np2
        if (idim>2) then
          j = jmax+np0
          do i=imin,imax,2
            data(i) = data(j)
            data(i+1) = -data(j+1)
            j = j-2
          end do
        endif
        j = jmax
        do i=imin,imax,np0
          data(i) = data(j)
          data(i+1) = -data(j+1)
          j = j-np0
        end do
      end do
    end do
!---
!--- END OF LOOP ON EACH DIMENSION
!---
900 np0 = np1
    np1 = np2
    nprev = n
  end do dimloop
end subroutine fft_gpd

!================================================
subroutine init_fft2_plan(plan)
  ! sets initial values of the plan structure
  type(planck_fft2_plan), intent(inout) :: plan

  plan%direction = fft2_forward
  plan%length    = -1
end subroutine init_fft2_plan


subroutine sanity_check (plan, len)
  type(planck_fft2_plan), intent(in) :: plan
  integer, intent(in) :: len

  if (len/=plan%length) &
    call exit_with_status(1,"planck_fft: invalid plan for this transform")
end subroutine sanity_check

subroutine make_fft2_plan (plan,length,direction)
  type(planck_fft2_plan), intent(out) :: plan
  integer, intent(in) :: length
  logical, intent(in) :: direction

  plan%length=length
  plan%direction=direction
end subroutine make_fft2_plan

subroutine destroy_fft2_plan (plan)
  type(planck_fft2_plan), intent(out) :: plan

  plan%direction=fft2_forward
  plan%length=-1
end subroutine destroy_fft2_plan

subroutine d_c_complex_fft2 (plan, data)
  ! replaced TRANSFER functions by loops for compatibility with GFORTRAN
  type(planck_fft2_plan), intent(in) :: plan
  complex(dp), intent(inout) :: data(:)

  real(dp) data2(2*size(data))
  integer :: i, lb

  lb = lbound(data,   dim = 1)
  call sanity_check (plan, size(data))
  do i=1,size(data)
    data2(2*i-1)=real (data(i+lb-1),kind=dp)
    data2(2*i)  =aimag(data(i+lb-1))
  end do
  call fft_gpd (data2,(/size(data)/),plan%direction,.false.)
  do i=1,size(data)
    data(i+lb-1) = cmplx(data2(2*i-1), data2(2*i), kind=dp)
  end do
end subroutine d_c_complex_fft2

subroutine d_r_complex_fft2 (plan, data)
  type(planck_fft2_plan), intent(in) :: plan
  real(dp), intent(inout) :: data(:)

  call sanity_check (plan, size(data)/2)
  call fft_gpd (data,(/size(data)/2/),plan%direction,.false.)
end subroutine d_r_complex_fft2

subroutine d_real_fft2 (plan, data)
  type(planck_fft2_plan), intent(in) :: plan
  real(dp), intent(inout) :: data(:)

  real(dp) data2(2*size(data))
  integer n, i

  call sanity_check (plan, size(data))
  n = size(data)

  if (plan%direction .eqv. fft2_forward) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),plan%direction,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),plan%direction,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine d_real_fft2

subroutine s_real_fft2 (plan, data)
  type(planck_fft2_plan), intent(in) :: plan
  real(sp), intent(inout) :: data(:)

  real(dp) data2(2*size(data))
  integer n, i

  call sanity_check (plan, size(data))
  n = size(data)

  if (plan%direction .eqv. fft2_forward) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),plan%direction,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),plan%direction,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine s_real_fft2
!======================================================================

subroutine complex_fft_orig (data, backward, onlyreal)
  complex(dp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward, onlyreal

  real(dp) data2(2*size(data))
  logical or, bw
  integer :: i, lb

  or = .false.
  if (present(onlyreal)) or=onlyreal
  bw = .false.
  if (present(backward)) bw=backward

  lb = lbound(data,   dim = 1)
!  data2 = transfer (data, data2, size=size(data2))
  do i=1,size(data)
     data2(2*i-1) = real (data(i+lb-1), kind=dp)
     data2(2*i  ) = aimag(data(i+lb-1))
  enddo
  call fft_gpd (data2,(/size(data)/),bw,or)
!  data = transfer (data2, data, size=size(data))
  do i=1,size(data)
     data(i+lb-1) = cmplx(data2(2*i-1), data2(2*i), kind=dp)
  enddo
end subroutine complex_fft_orig

subroutine complex_fft_alt (data, backward, onlyreal)
  real(dp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward, onlyreal

  logical or, bw

  or = .false.
  if (present(onlyreal)) or=onlyreal
  bw = .false.
  if (present(backward)) bw=backward

  call fft_gpd (data,(/size(data)/2/),bw,or)
end subroutine complex_fft_alt

subroutine d_real_fft (data, backward)
  real(dp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward

  real(dp) data2(2*size(data))
  logical bw
  integer n, i

  n = size(data)
  bw = .false.
  if (present(backward)) bw=backward

  if (.not. bw) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),bw,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),bw,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine d_real_fft

subroutine s_real_fft (data, backward)
  real(sp), intent(inout) :: data(:)
  logical, intent(in), optional :: backward

  real(dp) data2(2*size(data))
  logical bw
  integer n, i

  n = size(data)
  bw = .false.
  if (present(backward)) bw=backward

  if (.not. bw) then
    data2 = 0
    data2 (1:2*n-1:2) = data
    call fft_gpd (data2,(/n/),bw,.true.)
    data(1) = data2(1)
    data(2:n) = data2(3:n+1)
  else
    data2 = 0
    data2(1) = data(1)
    data2(3:n+1) = data(2:n)
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call fft_gpd (data2, (/n/),bw,.false.)
    data = data2(1:2*n-1:2)
  endif
end subroutine s_real_fft
!==================================================================
end module healpix_fft


module alm_tools
  !   Scalar+Open_MP implementation
  !
  !   function do_opemp
  !   subroutine init_rescale
  !   subroutine get_pixel_layout
  !   subroutine select_rings
  !   subroutine gen_recfac
  !   subroutine gen_recfac_spin
  !   subroutine gen_lamfac
  !   subroutine gen_lamfac_der
  !   subroutine gen_mfac
  !   subroutine gen_mfac_spin
  !   subroutine compute_lam_mm
  !   subroutine do_lam_lm
  !   subroutine do_lam_lm_spin
  !   subroutine do_lam_lm_pol
  !   subroutine gen_normpol
  !   function l_min_ylm
  !
  !   subroutine ring_synthesis
  !   subroutine ring_analysis
  !
  !   -------------------- in include files (see alm_map_dd_inc.f90) ---------
  !   subroutine alm2map_spin
  !   subroutine alm2map_pol
  !   subroutine alm2map_pol_pre1
  !   subroutine alm2map_pol_pre2
  !   subroutine alm2map_pol_der
  !
  !   subroutine map2alm_spin
  !   subroutine map2alm_pol
  !   subroutine map2alm_pol_pre1
  !   subroutine map2alm_pol_pre2
  !
  !   subroutine rotate_alm
  !   ------------------------------------------------------------------
  !
  !   subroutine plm_gen
  !
  use misc_utils, only: assert_alloc, assert, fatal_error
  use healpix_types
  use healpix_fft, only: real_fft2, planck_fft2_plan, make_fft2_plan, destroy_fft2_plan, fft2_forward, fft2_backward
  IMPLICIT none

  ! keep everything private unless stated otherwise
  private
  !--------------------------------------------------------------------------
  ! define large and small numbers used to renormalise the recursion on the Legendre Polynomials
  integer(I4B),      private, parameter :: LOG2LG   = 100
  real(KIND=DP),     private, parameter :: FL_LARGE = 2.0_dp **   LOG2LG
  real(KIND=DP),     private, parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
  ! declare array for dynamic rescaling of the Ylm
  integer(kind=i4b), private, parameter :: RSMAX = 20, RSMIN = -20
  real(dp),          private, dimension(RSMIN:RSMAX) :: rescale_tab
  real(DP),          private, parameter :: ALN2_INV = 1.4426950408889634073599246810_dp ! 1/log(2)
  ! misc
  integer(i4b),      private, parameter :: SMAXCHK = 50 ! maximum size of chunk (in number of ring pairs)
  ! parameters of Ylm short-cut
  integer(kind=i4b), private, parameter :: HPX_MXL0 = 40 ! minimum used, choose <=0 to do full loop
  real   (kind=dp),  private, parameter :: HPX_MXL1 = 1.35_dp
  !--------------------------------------------------------------------------

  ! make (front end) routines public
  public :: alm2map, map2alm, alm2map_der, alm2map_spin, map2alm_spin
  public :: rotate_alm
  public :: plm_gen
  public :: ring_synthesis, ring_analysis

  interface rotate_alm
     module procedure rotate_alm_d
  end interface
 
  interface alm2map_der
     module procedure alm2map_sc_der_d, alm2map_pol_der_d
  end interface

  interface alm2map_spin
     module procedure alm2map_spin_d
  end interface

  interface map2alm_spin
     module procedure map2alm_spin_d
  end interface

  interface alm2map
     module procedure alm2map_sc_d, alm2map_sc_pre_d, alm2map_pol_d, alm2map_pol_pre1_d, alm2map_pol_pre2_d
  end interface

  interface map2alm
     module procedure map2alm_sc_d, map2alm_sc_pre_d, map2alm_pol_d, map2alm_pol_pre1_d, map2alm_pol_pre2_d
  end interface

  interface sub_map2ring
     module procedure sub_map2ring_1d_s, sub_map2ring_1d_d, sub_map2ring_nd_s, sub_map2ring_nd_d
  end interface

  interface sub_ring2map
     module procedure sub_ring2map_1d_s, sub_ring2map_1d_d, sub_ring2map_nd_s, sub_ring2map_nd_d
  end interface

  ! make routines public as most of them are called by mpi_alm* routines
  public :: do_openmp
  public :: init_rescale
  public :: do_lam_lm, do_lam_lm_pol
  ! needed by mpi_alm*
  public :: l_min_ylm
  public :: get_pixel_layout
  public :: gen_recfac, gen_lamfac, gen_lamfac_der, gen_mfac, compute_lam_mm, gen_normpol
  public :: select_rings ! added Feb 2006
  !
  public :: gen_recfac_spin, gen_mfac_spin, do_lam_lm_spin ! July 2007

contains

  !**************************************************************************
  !             ALM2MAP/MAP2ALM    SIDEKICKS
  !**************************************************************************
  
function do_openmp()
    !================================================
    ! returns .true. if code was compiled with OpenMP
    !================================================
    logical(LGT) :: do_openmp
    !------------------------------------------------

    do_openmp = .false.
! DO NOT REMOVE FOLLOWING LINES
!$    do_openmp = .true.  ! Intel f90
!IBMP do_openmp = .true.  ! IBM xlf90
! -----------------------------

    return
end function do_openmp
  
  !================================================
  subroutine init_rescale()
    !================================================
    ! local variables
    integer(i4b) :: s, smax
    real(dp) :: logOVFLOW
    character(len=*), parameter :: code = 'gen_rescale'
    !------------------------------------------------
    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    if (smax > (RSMAX-1)) then
       print*,'Array rescale_tab too small in '//code
       print*,smax ,'>', RSMAX
       stop
    endif

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    do s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    enddo
    rescale_tab(0) = 1.0_dp

    return
  end subroutine init_rescale

  !=======================================================================
  subroutine get_pixel_layout(nside, ith, cth, sth, nphi, startpix, kphi0)
    !=======================================================================
    ! output Healpix pixel layout for the ring ith in [0,2*nside]
    !=======================================================================
    integer(I4B), intent(IN)  :: nside, ith
    real(DP)    , intent(OUT) :: cth, sth
    integer(I4B), intent(OUT) :: nphi, kphi0
    integer(I8B), intent(OUT) :: startpix
    !
    integer(I4B) :: nrings
    real(DP)     :: dth1, dth2, dst1
    !=======================================================================

    nrings = 2*nside
    if (ith < 1 .or. ith> nrings) then
       print*,'ith out of bounds ',ith,1,nrings
       call fatal_error
    endif

    dth1 = 1.0_dp / (3.0_dp*DBLE(nside)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nside))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nside) )

    if (ith < nside) then  ! polar cap (north)
       cth = 1.0_dp  - DBLE(ith)**2 * dth1
       nphi = 4*ith
       kphi0 = 1
       sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       startpix = 2*ith*(ith-1_I8B)
    else                   ! tropical band (north) + equator
       cth = DBLE(2*nside-ith) * dth2
       nphi = 4*nside
       kphi0 = MOD(ith+1-nside,2)
       sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       startpix = 2*nside*(nside-1_I8B) + (ith-nside)*int(nphi,kind=I8B)
    endif

    return
  end subroutine get_pixel_layout
  !=======================================================================
  subroutine select_rings(z, zbounds, keep_north, keep_south, keep_either)
    !=======================================================================
    ! select rings laying within zbounds
    ! if zbounds(1) < zbounds(2) : keep  zbounds(1) < z < zbounds(2)
    ! if zbounds(2) < zbounds(1) : keep z < zbounds(2) Union  zbounds(1) < z
    ! if zbounds(1)=zbounds(2) : keep everything
    ! input z should be >= 0
    !=======================================================================
    real(DP)    , intent(in)  :: z
    real(DP)    , intent(in), dimension(1:2)  :: zbounds
    logical(LGT), intent(OUT) :: keep_north, keep_south, keep_either
    !
    real(DP) :: zn, zs
    !=======================================================================

    ! if (zbounds(1) = zbounds(2)) keep everything
    if (abs(zbounds(1)-zbounds(2)) < 1.e-6) then
       keep_north    = .true.
       keep_south    = .true.
       keep_either   = .true.
       return
    endif

    zn = abs(z)
    zs = -zn

    if (zbounds(1) < zbounds(2)) then
       ! inner ring
       keep_north = (zn >= zbounds(1) .and. zn <= zbounds(2))
       keep_south = (zs >= zbounds(1) .and. zs <= zbounds(2))

    else
       ! outter ring
       keep_north = (zn > zbounds(1) .or. zn < zbounds(2))
       keep_south = (zs > zbounds(1) .or. zs < zbounds(2))
    endif
    keep_either   = keep_north .or. keep_south


    return
  end subroutine select_rings

  !=======================================================================
  subroutine gen_recfac( l_max, m, recfac)
  !=======================================================================
    ! generates recursion factors used to computes the Ylm of degree m 
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                            :: l_max, m
    real(DP),     intent(OUT), dimension(0:1, 0:l_max)  :: recfac
    !
    real(DP) :: fm2, fl2
    integer(I4B) :: l

    recfac(0:1,0:m) = 0.0_dp
    fm2 = DBLE(m) **2
    do l = m, l_max
       fl2 = DBLE(l+1) **2
       recfac(0,l) = DSQRT( (4.0_dp * fl2 - 1.0_dp) / (fl2-fm2) )
    enddo
    ! put outside the loop because of problem on some compilers
    recfac(1,m:l_max) = 1.0_dp / recfac(0,m:l_max)


    return
  end subroutine gen_recfac

  !=======================================================================
  subroutine gen_recfac_spin( l_max, m, spin, recfac_spin)
  !=======================================================================
    ! generates recursion factors used to computes the Ylms of degree m and spin s
    ! for all l in max(m,s)<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                            :: l_max, m, spin
    real(DP),     intent(OUT), dimension(0:2, 0:l_max)  :: recfac_spin
    !
    real(DP) :: fm2, fl2, fs2
    integer(I4B) :: l, l_low, aspin


    aspin = abs(spin)
    l_low = max(m, aspin)

    recfac_spin(0:2, 0:l_max) = 0.0_dp
    fm2 = DBLE(m) **2
    fs2 = DBLE(spin) **2
    do l = l_low, l_max
       fl2 = DBLE(l+1) **2
       recfac_spin(0,l) = DSQRT( (4.0_dp * fl2 - 1.0_dp) / (fl2-fm2) /(1.0_dp-fs2/fl2) )
    enddo
    do l = max(l_low, 1), l_max
       recfac_spin(2,l) = aspin * m / dble( l * (l+1) )
    enddo
    ! put outside the loop because of problem on some compilers
    recfac_spin(1,l_low:l_max) = 1.0_dp / recfac_spin(0,l_low:l_max)


    return
  end subroutine gen_recfac_spin

  !=======================================================================
  subroutine gen_lamfac( l_max, m, lam_fact)
  !=======================================================================
    ! generates factor relating scalar Ylm to polar Ylm
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: l_max, m
    real(DP),     intent(OUT), dimension(0:l_max)  :: lam_fact
    !
    real(DP) :: fm2, fl, fl2
    integer(I4B) :: l

!    lam_fact(0:m) = 0.
    lam_fact(0:max(1,m)) = 0.0_dp
    fm2 = real(m * m, kind=DP)
    do l = max(2,m+1), l_max
       fl  = real(l, kind=dp)
       fl2 = fl * fl
       lam_fact(l) = 2.0_dp * SQRT( (2.0_dp * fl + 1.0_dp) / (2.0_dp * fl - 1.0_dp) * (fl2-fm2) )
    enddo
    
    return
  end subroutine gen_lamfac

  !=======================================================================
  subroutine gen_lamfac_der(l_max, m, lam_fact)
    !=======================================================================
    ! generates factor relating scalar Ylm to its derivatives
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: l_max, m
    real(DP),     intent(OUT), dimension(0:l_max)  :: lam_fact
    !
    real(DP) :: fm2, fl, fl2
    integer(I4B) :: l

    lam_fact(0:m) = 0.
    fm2 = real(m * m, kind=DP)
    do l = max(1,m+1), l_max ! different lower bound than pol. factor
       fl  = real(l, kind=dp)
       fl2 = fl * fl
       lam_fact(l) = SQRT( (2.0_dp * fl + 1.0_dp) / (2.0_dp * fl - 1.0_dp) * (fl2-fm2) )
       ! different normalization than polarization factor
    enddo
    
    return
  end subroutine gen_lamfac_der

  !=======================================================================
  subroutine gen_mfac( m_max, m_fact)
  !=======================================================================
    ! generates factor used in lam_mm calculation
    ! for all m in 0<=m<=m_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: m_max
    real(DP),     intent(OUT), dimension(0:m_max)  :: m_fact
    !
    integer(I4B) :: m

    ! fact(m) = fact(m-1) * sqrt( (2m+1) / (2m) )
    m_fact(0) = 1.0_dp
    do m=1,m_max
      m_fact(m) = m_fact(m-1)*sqrt(dble(2*m+1)/dble(2*m))
    end do

    ! Log_2 ( fact(m) / sqrt(4 Pi) )
    do m=0,m_max
       m_fact(m) = log(SQ4PI_INV * m_fact(m)) * ALN2_INV 
    enddo

    return
  end subroutine gen_mfac
  !=======================================================================
  subroutine gen_mfac_spin( m_max, spin, m_fact)
  !=======================================================================
    ! generates factor used in lam_mm,s calculation
    ! for all m in 0<=m<=m_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: m_max, spin
    real(DP),     intent(OUT), dimension(0:m_max)  :: m_fact
    !
    integer(I4B) :: m, aspin
    real(DP) :: tmp

    aspin = abs(spin)
    m_fact(:) = -1.e30
    if (aspin <= m_max) m_fact(aspin) = 1.0_dp
    ! fact(m) = fact(m+1) * sqrt((s+m+1)/(s-m) )
    if (aspin > 0) then
       tmp = 1.d0
       do m = aspin-1, 0, -1
          tmp = tmp * sqrt( dble(aspin+m+1)/dble(aspin-m) )
          if (m <= m_max) m_fact(m) = tmp
       enddo
    endif
    ! fact(m) = fact(m-1) * sqrt(m*(2m+1)/(m-s)/(m+s)/2 )
    do m = aspin+1, m_max
      m_fact(m) = m_fact(m-1)*sqrt( m*dble(2*m+1)/dble(2*(m-aspin)*(m+aspin)) )
    end do

    ! Log_2 ( fact(m)  / sqrt(4 Pi))
    do m=0,m_max
       m_fact(m) = log(SQ4PI_INV * m_fact(m)) * ALN2_INV 
    enddo

    return
  end subroutine gen_mfac_spin
  !=======================================================================
  subroutine compute_lam_mm(mfac, m, sth, lam_mm, corfac, scalem)
    !=======================================================================
    ! computes lam_mm
    ! the true lam_mm is     lam_mm * corfac
    !=======================================================================
    integer(I4B),            intent(in)  :: m
    real(DP),                intent(in)  :: sth, mfac
    real(DP),                intent(out) :: lam_mm, corfac
    integer(I4B),            intent(out) :: scalem
    !
    real(DP) :: log2val, dlog2lg

    dlog2lg = real(LOG2LG, kind=DP)

    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalem = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalem,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalem * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m

    return
  end subroutine compute_lam_mm
  
  !=======================================================================
  subroutine do_lam_lm(lmax, m, cth, sth, mfac, recfac, lam_lm)
    !=======================================================================
    ! computes scalar lambda_lm(theta) for all l in [m,lmax] for a given m, and given theta
    ! input: lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        recfac: precomputed (by gen_recfac) quantities useful for 
    !            lambda_lm recursion for a given m
    ! output: lam_lm
    ! the routine also needs the array rescale_tac initialized by init_rescale
    !=======================================================================
    integer(I4B),                    intent(in)  :: lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(    0:lmax), intent(out) :: lam_lm
    !
    real(DP) :: log2val, dlog2lg
    real(DP) :: ovflow, unflow, corfac
    real(DP) :: lam_mm, lam_0, lam_1, lam_2
    integer(I4B) :: scalel, l, l_min
    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    l_min = l_min_ylm(m, sth)
    dlog2lg = real(LOG2LG, kind=DP)
    
    ! computes lamba_mm
    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalel = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalel,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalel * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m
    
    lam_lm(0:lmax) = 0.0_dp
    ! --- l = m ---
    lam_lm(m) = lam_mm * corfac

    ! --- l > m ---
    lam_0 = 0.0_dp
    lam_1 = 1.0_dp
    lam_2 = cth * lam_1 * recfac(0,m)
    do l = m+1, lmax
       ! do recursion
       if (l >= l_min) then
          lam_lm(l) = lam_2 * corfac * lam_mm
       endif
       lam_0 = lam_1 * recfac(1,l-1)
       lam_1 = lam_2
       lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

       ! do dynamic rescaling
       if (abs(lam_2) > ovflow) then
          lam_1 = lam_1*unflow
          lam_2 = lam_2*unflow
          scalel= scalel + 1
          corfac = rescale_tab(max(scalel,RSMIN))
       elseif (abs(lam_2) < unflow .and. abs(lam_2) /= 0.0_dp) then
          lam_1 = lam_1*ovflow
          lam_2 = lam_2*ovflow
          scalel= scalel - 1
          corfac = rescale_tab(max(scalel,RSMIN))
       endif
                   
    enddo ! loop on l
  end subroutine do_lam_lm
  !=======================================================================
  subroutine do_lam_lm_spin(lmax, m, spin, cth, sth, mfac, mfac_spin, recfac_spin, lam_lm_spin)
    !=======================================================================
    ! computes spin lambda_lm(theta) for all l in [m,lmax] for a given m, spin, and theta
    ! input: lmax, m, spin, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        mfac_spin: precomputed (by gen_mfac_spin) quantity useful for lambda_mm calculation
    !        recfac_spin: precomputed (by gen_recfac_spin) quantities useful for 
    !            lambda_lm recursion for a given m
    ! output: lam_lm
    ! the routine also needs the array rescale_tab initialized by init_rescale
    !=======================================================================
    integer(I4B),                    intent(in)  :: lmax,  m, spin
    real(DP),                        intent(in)  :: cth, sth, mfac, mfac_spin
    real(DP), dimension(0:2,0:lmax), intent(in)  :: recfac_spin
    real(DP), dimension(1:2,0:lmax), intent(out) :: lam_lm_spin
    !
    real(DP) :: dlog2lg
    real(DP) :: ovflow, unflow
    real(DP) :: lam_0, lam_1, lam_2
    real(DP) :: thetao2, ttho2, stho2, ctho2, kss, tmp1, tmp2
    real(DP), dimension(1:2) :: log10sss, corfac, log2val, lam_mm
    integer(I4B) :: l, l_min, l_low, aspin, iss
    integer(I4B), dimension(1:2) :: scalel
    !---------------------------------------------------------------

    lam_lm_spin(1:2,0:lmax) = 0.0_dp
    aspin = abs(spin)
    l_low = max(m, aspin)
    if (l_low > lmax) return

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    l_min = l_min_ylm(m, sth)
    dlog2lg = real(LOG2LG, kind=DP)
    
    thetao2 = atan2(sth, cth)/2.d0
    ttho2 = tan(thetao2)
    stho2 = sin(thetao2) ! sqrt(1-cth)/sqrt(2)
    ctho2 = cos(thetao2) ! sqrt(1+cth)/sqrt(2)
    ! lambda(s,s, s) = (-)^s * sqrt(2s+1/4Pi) sin(theta/2)^(2s)
    ! lambda(s,s,-s) = (-)^s * sqrt(2s+1/4Pi) cos(theta/2)^(2s)
    log10sss(1) = (2*aspin *log(stho2) + 0.5d0*log(2*aspin+1.0_dp)) ! log_10(abs(lam_sss)*sqrt(4Pi))
    log10sss(2) = (2*aspin *log(ctho2) + 0.5d0*log(2*aspin+1.0_dp)) ! log_10(abs(lam_ss-s)*sqrt(4Pi))
    !
    if (m >= aspin) then
       ! computes directly lambda(m,m,s)
       ! lambda(m,m,s) = - lambda(m-1,m-1,s) * sin(theta) * sqrt(m(2m+1)/2/(m^2-s^2))
       log2val(1:2) = mfac_spin + ((m-aspin)*log(sth) + log10sss(1:2)) * ALN2_INV ! log_2(lam_mms)
       scalel(1:2) = int (log2val(1:2) / dlog2lg)
       corfac(1) = rescale_tab(max(scalel(1),RSMIN))
       corfac(2) = rescale_tab(max(scalel(2),RSMIN))
       lam_mm(1:2) = 2.0_dp **(log2val(1:2) - scalel(1:2) * dlog2lg) ! rescaled lam_mm
       if (IAND(m,1)>0) lam_mm(1:2) = -lam_mm(1:2) ! negative for odd m
    else
       ! computes lambda(s,m,s)
       ! lambda(s,m,s)  = - lambda(s,m+1, s) / tan(theta/2) * sqrt((s+m+1)/(s-m))
       ! lambda(s,m,-s) = + lambda(s,m+1,-s) * tan(theta/2) * sqrt((s+m+1)/(s-m))
       log2val(1) = mfac_spin + ( (m-aspin)*log(ttho2) + log10sss(1)) * ALN2_INV ! log_2(lam_mms)
       log2val(2) = mfac_spin + (-(m-aspin)*log(ttho2) + log10sss(2)) * ALN2_INV ! log_2(lam_mms)
       scalel(1:2) = int (log2val(1:2) / dlog2lg)
       corfac(1) = rescale_tab(max(scalel(1),RSMIN))
       corfac(2) = rescale_tab(max(scalel(2),RSMIN))
       lam_mm(1:2) = 2.0_dp **(log2val(1:2) - scalel(1:2) * dlog2lg) ! rescaled lam_mm
       if (IAND(m,1)>0)     lam_mm(1) = -lam_mm(1) ! negative for odd m, and s>0
       if (IAND(aspin,1)>0) lam_mm(2) = -lam_mm(2) ! negative for odd s, and s<0
    endif

    ! --- l = m ---
    lam_lm_spin(1:2,l_low) = lam_mm(1:2) * corfac(1:2)

    do iss=1,2 ! +spin and then -spin
       kss = 1.0_dp
       if (iss == 2) kss = -1.0_dp
       ! --- l > m ---
       lam_0 = 0.0_dp
       lam_1 = 1.0_dp
       lam_2 = (cth + kss*recfac_spin(2,l_low)) * lam_1 * recfac_spin(0,l_low)
       do l = l_low+1, lmax
          ! do recursion
          if (l >= l_min) then
             lam_lm_spin(iss,l) = lam_2 * corfac(iss) * lam_mm(iss)
          endif
          lam_0 = lam_1 * recfac_spin(1,l-1)
          lam_1 = lam_2
          lam_2 = ((cth + kss*recfac_spin(2,l)) * lam_1 - lam_0) * recfac_spin(0,l)

          ! do dynamic rescaling
          if (abs(lam_2) > ovflow) then
             lam_1 = lam_1*unflow
             lam_2 = lam_2*unflow
             scalel(iss)= scalel(iss) + 1
             corfac(iss) = rescale_tab(max(scalel(iss),RSMIN))
          elseif (abs(lam_2) < unflow .and. abs(lam_2) /= 0.0_dp) then
             lam_1 = lam_1*ovflow
             lam_2 = lam_2*ovflow
             scalel(iss)= scalel(iss) - 1
             corfac(iss) = rescale_tab(max(scalel(iss),RSMIN))
          endif
       enddo ! loop on l
    enddo ! loop on spin sign (iss)

    ! compute functions with fixed parity
    ! beware: the first one is always the half sum
    !         the second one is always the half difference, independently of spin 
    do l=0, lmax
       tmp1 = lam_lm_spin(1,l) * 0.5_dp
       tmp2 = lam_lm_spin(2,l) * 0.5_dp
       lam_lm_spin(1,l) = tmp1 + tmp2
       lam_lm_spin(2,l) = tmp1 - tmp2
    enddo

  end subroutine do_lam_lm_spin
  !=======================================================================
  subroutine do_lam_lm_pol(lmax, m, cth, sth, mfac, recfac, lam_fact, lam_lm)
    !=======================================================================
    ! computes temperature&polar lambda_lm(p,theta) for all l in [m,lmax] for a given m, and given theta
    ! input: lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        recfac: precomputed (by gen_recfac) quantities useful for 
    !            lambda_lm recursion for a given m
    !        lam_fact: precomputed (by gen_lamfac) factor useful for polarised lambda recursion
    ! output: lam_lm for T and P
    ! the routine also needs the array rescale_tac initialized by init_rescale
    !=======================================================================
    integer(I4B),                    intent(in)  :: lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(    0:lmax), intent(in)  :: lam_fact
    real(DP), dimension(1:3,0:lmax), intent(out) :: lam_lm
    !
    real(DP) :: log2val, dlog2lg
    real(DP) :: ovflow, unflow, corfac
    real(DP) :: lam_mm, lam_0, lam_1, lam_2, lam_lm1m
    integer(I4B) :: scalel, l, l_min
    real(DP) :: normal_m, fm2, fl, flm1
    real(DP) :: two_cth, one_on_s2, fm_on_s2, two_on_s2, c_on_s2
    real(DP) :: a_w, a_x, b_w
    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    l_min = l_min_ylm(m, sth)
    dlog2lg = real(LOG2LG, kind=DP)
    
    fm2       = real(m * m, kind = DP)
    normal_m  = (2.0_dp * m) * (1 - m)
    two_cth   = 2.0_dp * cth
    one_on_s2 = 1.0_dp / (sth * sth)
    fm_on_s2  =      m * one_on_s2
    two_on_s2 = 2.0_dp * one_on_s2
    c_on_s2   = cth    * one_on_s2
    b_w       =  c_on_s2 
    

    ! computes lamba_mm
    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalel = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalel,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalel * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m
    
    lam_lm(1:3,m:lmax) = 0.0_dp
    ! --- l = m ---
    lam_lm(1,m) = corfac * lam_mm !Actual lam_mm 

    if (m >= l_min) then ! skip Ymm if too small
       lam_lm(2,m) =  (normal_m * lam_lm(1,m))  * ( 0.5_dp - one_on_s2 )
       lam_lm(3,m) =  (normal_m * lam_lm(1,m))  *            c_on_s2
    endif

    ! --- l > m ---
    lam_0 = 0.0_dp
    lam_1 = 1.0_dp
    lam_2 = cth * lam_1 * recfac(0,m)

    do l = m+1, lmax
       ! do recursion
       lam_lm1m = lam_lm(1,l-1) * lam_fact(l) ! must be incremented even if not used
       lam_lm(1,l) = lam_2 * corfac * lam_mm
       if (l >= l_min) then
          fl = real(l, kind = DP)
          flm1 = fl - 1.0_dp
          a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
          a_x =  two_cth * flm1
          lam_lm(2,l) =                b_w * lam_lm1m - a_w * lam_lm(1,l)
          lam_lm(3,l) = fm_on_s2 * (         lam_lm1m - a_x * lam_lm(1,l))
       endif
       lam_0 = lam_1 * recfac(1,l-1)
       lam_1 = lam_2
       lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

       ! do dynamic rescaling
       if (abs(lam_2) > ovflow) then
          lam_1 = lam_1*unflow
          lam_2 = lam_2*unflow
          scalel= scalel + 1
          corfac = rescale_tab(max(scalel,RSMIN))
       elseif (abs(lam_2) < unflow .and. abs(lam_2) /= 0.0_dp) then
          lam_1 = lam_1*ovflow
          lam_2 = lam_2*ovflow
          scalel= scalel - 1
          corfac = rescale_tab(max(scalel,RSMIN))
       endif
                   
    enddo ! loop on l
  end subroutine do_lam_lm_pol
  !=======================================================================
  subroutine gen_normpol(l_max, normal_l)
    !=======================================================================
    ! generates normalisaton factors for polarisation basis functions
    ! for all l 
    !=======================================================================
    integer(I4B), intent(IN)                       :: l_max
    real(DP),     intent(OUT), dimension(0:l_max)  :: normal_l
    !
    integer(I4B) :: l
    real(DP)     :: fl, xx

    normal_l(0:1) = 0.0_dp
    do l = 2, l_max
       fl = DBLE(l)
       xx = (fl+2.0_dp) * (fl+1.0_dp) * fl * (fl-1.0_dp)
       normal_l(l) = SQRT ( KvS / xx)
       ! either CMBFAST (KvS=1) or Kamionkowski et al. (KvS=2) definition
    enddo

    return
  end subroutine gen_normpol

  !================================================================
  function l_min_ylm(m, sth) result(lmin)
  !================================================================
    ! returns minimal order l at which to keep Ylm
    ! |Ylm| < eps * Y00 ==>
    ! m_cut(theta, l) = theta * l * e / 2 + | ln(eps)| + ln(l)/2
    ! if eps = 1.e-15 and l < 1.e4
    ! m_cut(theta, l) = theta * l * 1.35 + 40
    ! the choice of 1.35 (or larger) 
    ! also insures that the equatorial rings will have all their Ylm's computed
    ! default parameters are HPX_MXL0 = 40 and HPX_MXL1 = 1.35_DP
    !======================================================
    ! parameters of short-cut: defined in module header
    ! dummy variables
    integer(I4B)             :: lmin
    integer(I4B), intent(IN) :: m
    real(DP),     intent(IN) :: sth

    lmin = m ! default
    if (HPX_MXL0 > 0) lmin = max(lmin, int((m - HPX_MXL0)/(HPX_MXL1 * sth)))

    return
  end function l_min_ylm
  !=========================================================
  !**************************************************************************
  !
  !             FOURIER TRANSFORM ON RINGS
  !
  !**************************************************************************
  !=======================================================================
  subroutine ring_synthesis(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)
    !=======================================================================
    !     RING_SYNTHESIS
    !       called by alm2map
    !       calls     real_fft
    !
    !     dataout(j) = sum_m datain(m) * exp(i*m*phi(j))
    !     with phi(j) = j*2pi/nph + kphi0*pi/nph and kphi0 =0 or 1
    !
    !     as the set of frequencies {m} is larger than nph,
    !     we wrap frequencies within {0..nph-1}
    !     ie  m = k*nph + m' with m' in {0..nph-1}
    !     then
    !     noting bw(m') = exp(i*m'*phi0)
    !                   * sum_k (datain(k*nph+m') exp(i*k*pi*kphi0))
    !        with bw(nph-m') = CONJ(bw(m')) (if datain(-m) = CONJ(datain(m)))
    !     dataout(j) = sum_m' [ bw(m') exp (i*j*m'*2pi/nph) ]
    !                = Fourier Transform of bw
    !        is real
    !
    !         NB nph is not necessarily a power of 2
    !
    !=======================================================================


    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph, kphi0
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(IN)  :: datain
    REAL(DP),     DIMENSION(0:nph-1), INTENT(OUT) :: dataout

    INTEGER(I4B) :: iw,ksign,m,k,kshift
    COMPLEX(DPC), DIMENSION(0:nph-1) :: bw
    type(planck_fft2_plan) :: plan
    COMPLEX(DPC) :: dat
    real(DP)     :: arg
    !=======================================================================

    ksign = + 1
    kshift = (-1)**kphi0  ! either 1 or -1
    bw(0:nph-1) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

    !     all frequencies [-m,m] are wrapped in [0,nph-1]
    bw(0)=datain(0)
    do m  = 1, nmmax                        ! in -nmmax, nmmax
       iw = MODULO(m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + datain(m)*(kshift**k)  ! complex number
       iw = MODULO(-m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (-m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + CONJG(datain(m))*(kshift**k)  ! complex number
    enddo

    !     kshift**k = 1       for even turn numbers
    !               = 1 or -1 for odd  turn numbers : results from the shift in space

    !     applies the shift in position <-> phase factor in Fourier space
    dataout(0)=REAL(bw(0), kind=DP)
    do iw = 1, nph/2-1
       m = ksign*(iw)
       if(kphi0==1) then
          arg = m * PI / dble(nph)
          dat =bw(iw) * CMPLX( DCOS(arg), DSIN(arg), kind=DP)
       else
          dat =bw(iw)
       endif
       dataout(iw*2-1) = REAL(dat, kind=DP)
       dataout(iw*2  ) = AIMAG(dat)
    enddo
    iw=nph/2
    m = ksign*(iw)
    if(kphi0==1) then
       arg = m * PI / dble(nph)
       dat =bw(iw) * CMPLX( DCOS(arg), DSIN(arg), kind=DP)
    else
       dat =bw(iw)
    endif
    dataout(iw*2-1) = REAL(dat, kind=DP)

    call make_fft2_plan(plan,length=nph,direction=fft2_backward)
    call real_fft2 (plan, dataout)
    call destroy_fft2_plan (plan)
    RETURN

  END subroutine ring_synthesis

  !=======================================================================
  subroutine ring_analysis(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)
    !=======================================================================
    !     ring_analysis
    !       called by map2alm
    !       calls     real_fft
    !
    !     integrates (data * phi-dependence-of-Ylm) over phi
    !     --> function of m can be computed by FFT
    !     with  0<= m <= npoints/2 (: Nyquist)
    !     because the data is real the negative m are the conjugate of the
    !     positive ones
    !=======================================================================

    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph, kphi0
    REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(OUT) :: dataout

    INTEGER(I4B) :: i,m,im_max,ksign
!    REAL(DP) :: phi0
    REAL(DP), DIMENSION(0:nph-1) :: data
    type(planck_fft2_plan) :: plan
    real(DP)  :: arg

    !=======================================================================

    ksign = - 1
    data=0.
    data(0:nph-1)=datain(0:nph-1)

    call make_fft2_plan(plan,length=nph,direction=fft2_forward)
    call real_fft2(plan,data)
    call destroy_fft2_plan(plan)

    im_max = MIN(nph/2,nmmax)

    ! m = 0,  i=0
    dataout(0)=CMPLX(data(0),0.0_dp,kind=DP)

    ! 1 <= m <= im_max, --> i=1,2,3,..., im_max
    do i = 1, im_max*2-3, 2
       dataout((i+1)/2) = CMPLX( data(i), data(i+1),kind= DP)
    enddo

    if(im_max==nph/2) then
       dataout(im_max)= CMPLX( data(nph-1),0.0_dp,kind=DP) ! m = Nyquist
    else
       dataout(im_max)= CMPLX( data(2*im_max-1),data(2*im_max),kind=DP)
    endif

    if(im_max==nmmax) goto 1000 ! m_max <= Nyquist

    ! if (m > Nyquist)
    do i =  im_max+1,min(nph,nmmax)
       dataout(i) = conjg(dataout(2*im_max-i) )
    end do

    if(min(nph,nmmax)==nmmax) goto 1000 ! nph > nmmax

    do i =  2*im_max+1,nmmax
       dataout(i) = dataout(mod(i,2*im_max))
    end do

1000 continue

    if(kphi0==1)then
       do i =  0,nmmax
          m = ksign*i
!           dataout(i)=dataout(i)* CONJG(trig(-m,nph/4))
          arg = m * PI / dble(nph)
          dataout(i)=dataout(i)* CMPLX( DCOS(arg), DSIN(arg), kind=DP)
       enddo
    end if

    RETURN
  END subroutine ring_analysis

  !**************************************************************************
  !   ALM2MAP
  !   MAP2ALM
  !**************************************************************************

  ! single precision routines
!#include "alm_map_ss_inc.F90"

  ! double precision routines
#include "alm_map_dd_inc.F90"
  
  !**************************************************************************
  !             PLM GENERATION
  !**************************************************************************
  !========================================================
  subroutine plm_gen(nsmax, nlmax, nmmax, plm)
    !========================================================
    use long_intrinsic, only: long_size

    integer(i4b),             intent(IN) :: nsmax, nlmax, nmmax
    real(dp), dimension(0:,1:), intent(OUT):: plm

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nd2, nrings
    integer(I8B) :: nd1, n_lm, n_plm, i_mm, i_up
    real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
!!    real(DP)                 :: lambda_w, lambda_x
    real(DP)                 :: normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1, fpol
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)                 :: ovflow, unflow
    real(DP),     dimension(:,:,:), allocatable :: plm_sub
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP),     dimension(:),   allocatable :: lam_fact
    real(DP),     dimension(:), allocatable :: mfac

    real(DP),     dimension(:),     allocatable :: normal_l
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl
    integer(i4b)        :: nph, kphi0, i
    integer(i8b)        :: startpix
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2

    INTEGER(I4B) :: status
    LOGICAL(LGT) :: polarisation
    character(len=*), parameter :: code = 'PLM_GEN'

    !=================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat

    n_lm  = ((nmmax+1)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm = n_lm * nrings
    nd1 = long_size(plm, 1)
    nd2 = long_size(plm, 2)

    if (nd1 < n_plm) then
       print*,code//' > Plm array too small:', nd1, n_plm
       stop
    endif
    if (nd2 /= 1 .and. nd2 /= 3) then
       print*,code//' > Plm array should have dimension 1 or 3:',nd2
       stop
    endif
    polarisation = (nd2 == 3)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')
    if (polarisation) then
       allocate(normal_l(0:nlmax),stat = status)
       call assert_alloc(status,code,'normal_l')
    endif

    if (.not. do_openmp()) then
       allocate(recfac(0:1,0:nlmax), plm_sub(1:nd2,0:nlmax,0:chunksize-1), stat = status)
       call assert_alloc(status,code,'recfac & plm_sub')
       if (polarisation) then
          allocate(lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'lam_fact')
       endif
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax, mfac)
    ! generate Polarization normalisation
    if (polarisation) call gen_normpol(nlmax, normal_l)
    call init_rescale()
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    plm = 0.0_dp

    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph, startpix, kphi0)
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
       enddo

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, nd2, chunksize, n_lm, &
!$OMP    rescale_tab, ovflow, unflow, polarisation, &
!$OMP    plm, cth, sth, mfac, normal_l, one_on_s2, c_on_s2) &
!$OMP private(plm_sub, recfac, lam_fact, status, &
!$OMP   m, fm2, normal_m, ithl, ith, i_mm, i_up, &
!$OMP   scalem, scalel, corfac, fpol, &
!$OMP   lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2, &
!$OMP   cth_ring, fm_on_s2, two_on_s2, two_cth_ring, a_w, a_x, b_w,  &
!$OMP   l, fl, flm1, i)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(recfac(0:1,0:nlmax), plm_sub(1:nd2,0:nlmax,0:chunksize-1), stat = status)
          call assert_alloc(status,code,'recfac & plm_sub')
          if (polarisation) then
             allocate(lam_fact(0:nlmax),stat = status)
             call assert_alloc(status,code,'lam_fact')
          endif
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          if (polarisation) then
             ! generate Ylm relation factor for degree m
             call gen_lamfac(nlmax, m, lam_fact)
             fm2 = real(m * m, kind = DP)
             normal_m = (2.0_dp * m) * (1 - m)
          endif

          do ithl = 0, uchk - lchk
             ! determine lam_mm
             call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)
             ! ---------- l = m ----------
             !           temperature 
             lam_lm = corfac * lam_mm !Actual lam_mm 
             plm_sub(1, m, ithl) = lam_lm

             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             cth_ring = cth(ithl)
             lam_2 = cth_ring * lam_1 * recfac(0,m)

             if (polarisation) then
                fpol = normal_m * normal_l(m) * lam_lm
                plm_sub(2, m, ithl) =  fpol * ( 0.5_dp - one_on_s2(ithl) )
                plm_sub(3, m, ithl) =  fpol *              c_on_s2(ithl)
                !
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
             endif
             ! ---------- l > m ----------
             do l = m+1, nlmax
                if (polarisation) lam_lm1m = lam_lm * lam_fact(l)
                lam_lm = lam_2 * corfac * lam_mm
                plm_sub(1, l, ithl) = lam_lm
                
                if (polarisation) then
                   fl = real(l, kind = DP)
                   flm1 = fl - 1.0_dp
                   a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                   a_x =  two_cth_ring * flm1
                   plm_sub(2, l, ithl) =            (   b_w * lam_lm1m - a_w * lam_lm) * normal_l(l)
                   plm_sub(3, l, ithl) = fm_on_s2 * (         lam_lm1m - a_x * lam_lm) * normal_l(l)
                endif

                lam_0 = lam_1 * recfac(1,l-1)
                lam_1 = lam_2
                lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                if (abs(lam_2) > OVFLOW) then
                   lam_1 = lam_1*UNFLOW
                   lam_2 = lam_2*UNFLOW
                   scalel= scalel + 1
                   corfac = rescale_tab(max(scalel+scalem,RSMIN))
                elseif (abs(lam_2) < UNFLOW) then
                   lam_1 = lam_1*OVFLOW
                   lam_2 = lam_2*OVFLOW
                   scalel= scalel - 1
                   corfac = rescale_tab(max(scalel+scalem,RSMIN))
                endif
             enddo ! loop on l
          enddo ! loop on rings (ithl)

          ! do memory skipping operations outside inner loops
          do ith = lchk, uchk
             i_mm = n_lm * (ith-1) + ((2*nlmax+3-m)*m)/2 ! location of Ym,m for ring ith
             i_up = i_mm + nlmax - m ! location of Ynlmax,m for ring ith
             ithl = ith - lchk
             do i=1, nd2
                plm(i_mm:i_up, i) = plm_sub(i, m:nlmax, ithl)
             enddo
          enddo
!!!!!!!!
!           i_mm = n_lm*chunksize*ichunk + ((2*nlmax+3-m)*m)/2
!           do ithl = 0, uchk-lchk
!              i_up = i_mm+(nlmax-m+1)*ithl
!           do i=1,nd2
!              plm(i_up:i_up+nlmax-m,i) = plm_sub(i, m:nlmax, ithl)
!           enddo
!           enddo
          !------------------------------------------------
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (recfac, plm_sub)
          if (polarisation) deallocate(lam_fact)
       endif
!$OMP end parallel


    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate (recfac, plm_sub)
       if (polarisation) deallocate(lam_fact)
    endif
    deallocate(mfac)
    if (polarisation) deallocate(normal_l)

    return
  end subroutine plm_gen


end module alm_tools



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
module misc_utils
!   subroutine fatal_error
!   function file_present
!   subroutine assert_directory_present
!   subroutine assert_not_present
!   subroutine assert_alloc
!   subroutine assert
!   function expand_env_var
!   function string
!   subroutine brag_openmp
!--------------------------------------------------------------------------------------
  ! edited 2006-10-31: fatal_error for gfortran gcc4.1.1 bug workaround (V. Stolyarov)
  ! 2007-06-06: string() now accepts LGT variables
  ! 2008-03-25: added expand_env_var, string() accepts 64-bit integer variables (on systems that can deal with them)
  ! 2012-10-29: edited file_present to accept virtual files and CFITSIO 'extended filenames'
  !    the NOCFITSIO flag must be set to return to standard UNIX 'inquire' behavior
  !--------------------------------------------------------------------------------------
  use healpix_types
  use extension, only : exit_with_status, getEnvironment
  implicit none
  private

  integer, parameter, private :: LCH=48
  interface string
#ifdef NO64BITS
     module procedure string_l, string_i,           string_s, string_d
#else
     module procedure string_l, string_i, string_j, string_s, string_d
#endif
  end interface

  interface fatal_error
  	module procedure fatal_error_womsg, fatal_error_msg
  end interface


  public :: fatal_error, assert, assert_alloc
  public :: brag_openmp
  public :: expand_env_var

contains
  !-----------------------------------------------------
!  subroutine fatal_error (msg)
!    character(len=*), intent(in), optional :: msg
!
!    if (present(msg)) then
!       print *,'Fatal error: ', trim(msg)
!    else
!       print *,'Fatal error'
!    endif
!    call exit_with_status(1)
!  end subroutine fatal_error

  subroutine fatal_error_msg (msg)
    character(len=*), intent(in) :: msg
       print *,'Fatal error: ', trim(msg)
    call exit_with_status(1)
  end subroutine fatal_error_msg

  subroutine fatal_error_womsg
      print *,'Fatal error'
    call exit_with_status(1)
  end subroutine fatal_error_womsg

  !-----------------------------------------------------
  subroutine assert_alloc (stat,code,arr)
!!    integer, intent(in) :: stat
    integer(i4b), intent(in) :: stat
    character(len=*), intent(in) :: code, arr

    if (stat==0) return

    print *, trim(code)//'> cannot allocate memory for array: '//trim(arr)
    call exit_with_status(1)
  end subroutine assert_alloc

  !-----------------------------------------------------
  subroutine assert (testval,msg,errcode)
    logical, intent(in) :: testval
    character(len=*), intent(in), optional :: msg
    integer(i4b), intent(in), optional :: errcode

    if (testval) return

    print *,"Assertion failed: "
    if (present(msg)) print *, trim(msg)
    if (present(errcode)) call exit_with_status (errcode)
    call exit_with_status(1)
  end subroutine assert

  !-----------------------------------------------------
  function expand_env_var(instr) result(outstr)
    ! substitute an evironment variable invocation ${VAR} 
    ! with its value in a string
    character(len=*), intent(in) :: instr
    character(len=FILENAMELEN)   :: outstr, tmp, varname, varvalue
    integer  :: i1, i2, ln
    character(len=*), parameter :: code = 'expand_env_var'

    outstr = trim(adjustl(instr))

    ! replace leading ~/ with the value of $HOME
    i1 = index(outstr,'~/')
    if (i1 == 1) then
       ln = len_trim(outstr)
       call getEnvironment('HOME',varvalue)
       tmp = trim(adjustl(varvalue))//outstr(2:ln)
       outstr = trim(adjustl(tmp))
    endif

    ! replace ${XXX} with the value of $XXX
    do
       ln = len_trim(outstr)
       i1 = index(outstr,'${')
       if (i1 <= 0) exit
       i2 = index(outstr,'}')
       
       if (i2 <=  i1 + 1) then
          print*,'WARNING: '//code//' can not process string: '//trim(instr)
          print*,'         Unmatched  { or } .'
          outstr = instr
          return
       endif
       varname = outstr(i1+2:i2-1)
       call getEnvironment(varname, varvalue)
       tmp = outstr(1:i1-1)//trim(adjustl(varvalue))//outstr(i2+1:ln)
       outstr = trim(adjustl(tmp))
    enddo
    
    return
  end function expand_env_var

  subroutine brag_openmp()
  !================================================
    ! OpenMP bragging
    !================================================
    ! OpenMP variables
    !$     integer :: omp_get_thread_num, omp_get_num_threads, omp_get_num_procs
    !IBMP  integer :: omp_get_thread_num, omp_get_num_threads, omp_get_num_procs

!$OMP parallel
!
!$   if (omp_get_thread_num() == 0) then
!$       write(*,9000) ' --------------------------------------'
!$       write(*,9010) ' Number of OpenMP threads in use: ', omp_get_num_threads()
!$       write(*,9010) ' Number of CPUs available:        ', omp_get_num_procs()
!$       write(*,9000) ' --------------------------------------'
!$   end if
!
!IBMP   if (omp_get_thread_num() == 0) then
!IBMP       write(*,9000) ' --------------------------------------'
!IBMP       write(*,9010) ' Number of OpenMP threads in use: ', omp_get_num_threads()
!IBMP       write(*,9010) ' Number of CPUs available:        ', omp_get_num_procs()
!IBMP       write(*,9000) ' --------------------------------------'
!IBMP   end if
!
!$OMP end parallel
9000 format(a)
9010 format(a,i4)

    return
  end subroutine brag_openmp

end module misc_utils


module extension
  ! defines in F90 some C commands
  ! These extensions are not completely standard in all F90/95 compilers
  !
  ! getEnvironment   : emulates getenv
  ! getArgument      : emulates getarg
  ! nArguments       : emulates iargc
  !
  ! written by Eric Hivon, Nov 2001
  !
  ! exit_with_status : verbose and clean exit, added by M.R.
  ! 2005-08: edited for Gfortran
  ! 2013-05-07: G95-compatible
  ! 2015-07-31: G95-compatible
  ! 2016-05: edited for __GFORTRAN__
  
#ifdef NAG
  USE f90_unix, ONLY : iargc, getarg, exit
#endif
!VF  USE dflib, ONLY : nargs, getarg

  USE healpix_types, ONLY : I4B, I8B
  IMPLICIT none

#if ((!defined(NAG)) && (!defined(GFORTRAN)) && (!defined(__GFORTRAN__)))
interface
  function iargc()
    integer iargc
  end function

  subroutine getarg (num, res)
    integer, intent(in) :: num
    character(len=*), intent(out) :: res
  end subroutine
end interface
#endif

! work-around G95 bug (2013-05-07)
  integer(kind=I4B), parameter, private :: arg_shift = 0
!  integer(kind=I4B), private :: arg_shift = 0
!VF  integer(kind=I4B), private :: arg_shift = 1

#ifdef NO64BITS
  interface exit_with_status
     module procedure exit_with_status
  end interface
#else
  interface exit_with_status
     module procedure exit_with_status, exit_with_status_8
  end interface
#endif
  private
  public :: getEnvironment, getArgument, nArguments, exit_with_status

  contains

#if (defined (GFORTRAN) || defined(__GFORTRAN__) )

    ! ===========================================================
    function iargc ()
    ! ===========================================================
       integer iargc
       ! ===========================================================

       iargc=command_argument_count()
    end function

    ! ===========================================================
    subroutine getarg(num, res)
    ! ===========================================================
       !integer, intent(in) :: num ! G95, 2015-07-30
       integer(i4b), intent(in) :: num
       character(len=*), intent(out) :: res
       integer num1, l, err
       ! ===========================================================
       num1 = num
       call get_command_argument(num1,res,l,err)
    end subroutine

#endif

    ! ===========================================================
    function nArguments() result(narg)
      ! ===========================================================
      integer(kind=I4B) :: narg
      ! ===========================================================

      narg = iargc() - arg_shift
!VF      narg = nargs() - arg_shift

      return
    end function nArguments
    ! ===========================================================
    subroutine getEnvironment(name, value)
      ! ===========================================================
      character(len=*), intent(in) :: name
      character(len=*), intent(out) :: value
      integer(kind=I4B) :: inull, lnstr
!       character(len=200) :: name1
      ! ===========================================================
      ! call C routine after adding a trailing NULL to input
      value = ""
      call cgetenvironment(trim(adjustl(name))//char(0), value)
      ! remove trailing NULL (\0) created by C routine on output
      lnstr = len(value)
      inull = index(value, char(0), back=.true.)
      if (inull > 0 .and. inull < lnstr) value(inull:inull) = " "

      return
    end subroutine getEnvironment
    ! ===========================================================
    subroutine getArgument(index, argument)
      ! ===========================================================
      integer(kind=I4B), intent(in) :: index
      character(len=*), intent(out) :: argument
      integer(kind=I4B) :: i1
      ! ===========================================================
      i1 = index + arg_shift
      call getarg(i1, argument)

      return
    end subroutine getArgument

! i4b and i8b versions of exit_with_status   ! G95 2015-07-30
    ! ===========================================================
    subroutine exit_with_status (code, msg)
      ! ===========================================================
      integer(i4b), intent(in) :: code
      character (len=*), intent(in), optional :: msg
      ! ===========================================================
      if (present(msg)) print *,trim(msg)
      print *,'program exits with exit code ', code
#if (defined (RS6000))
      call exit_ (code)
#else
      call exit (code)
#endif
    end subroutine exit_with_status

#ifndef NO64BITS
    ! ===========================================================
    subroutine exit_with_status_8 (code, msg)
      ! ===========================================================
      integer(i8b), intent(in) :: code
      character (len=*), intent(in), optional :: msg
      ! ===========================================================
      if (present(msg)) print *,trim(msg)
      print *,'program exits with exit code ', code
#if (defined (RS6000))
      call exit_ (code)
#else
      call exit (code)
#endif
    end subroutine exit_with_status_8
#endif

end module extension

