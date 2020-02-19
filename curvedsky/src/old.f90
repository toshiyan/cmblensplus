!////////////////////////////////////////////////////!
! Map filtering
!////////////////////////////////////////////////////!

module cninv
  !from F90/src_utils
  use constants, only: pi
  use general, only: savetxt
  use hp_cgd
  implicit none

  private pi

contains


subroutine cnfilter(npix,lmax,cl,nij,alm,xlm,chn,lmaxs,nsides,itns,eps,fratio,filter,verbose)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*    1) The signal power spectrum is isotropic Gaussian. 
!*    2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*    3) The data model is bxS+N
!*
!* Args:
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:lmax)
!*    :nij[pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:npix-1)
!*    :alm[l,m] (dcmplx) : Input alm, with bouds (0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*
!* Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!* Returns:
!*    :xlm[l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: fratio
  integer, intent(in) :: npix, lmax, chn
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: fratio = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(0:lmax) :: cl
  double precision, intent(in), dimension(0:npix-1) :: nij
  double complex, intent(in), dimension(0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(0:lmax,0:lmax) :: xlm
  !internal
  type(mgchain) :: mgc
  integer :: ni, l, ilmaxs(chn), insides(chn)
  double precision :: clh(1:1,0:lmax,0:lmax), ninv(1:1,0:npix-1), ratio(itns(1))
  double complex :: b(1:1,0:lmax,0:lmax), nalm(1:1,0:lmax,0:lmax), nxlm(1:1,0:lmax,0:lmax)

  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides = (/int(dsqrt(npix/12d0))/)
  else
    ilmaxs  = lmaxs
    insides = nsides
  end if

  clh = 0d0
  do l = 1, lmax
    clh(1,l,0:l) = dsqrt(cl(l))
  end do

  !first compute b = C^1/2 N^-1 X
  ninv(1,:)   = nij
  nalm(1,:,:) = alm
  call mat_multi(1,npix,lmax,clh,ninv,nalm,b,'rhs')

  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call set_mgchain(mgc,chn,ilmaxs,insides,itns,eps,verbose)
  call cg_algorithm(1,clh,ninv,b,nxlm,mgc,1,ratio)
  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  !convert
  xlm = 0d0
  do l = 1, lmax
      if (filter=='')   xlm(l,0:l) = nxlm(1,l,0:l)/clh(1,l,0:l)
      if (filter=='W')  xlm(l,0:l) = nxlm(1,l,0:l)*clh(1,l,0:l)
  end do

end subroutine cnfilter


subroutine cnfiltertp(n,npix,lmax,cl,nij,alm,xlm,chn,lmaxs,nsides,itns,eps,fratio,filter,verbose)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is bxS+N
!*
!* Args:
!*    :n (int) : Number of maps
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :nij[n,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :alm[n,l,m] (dcmplx) : Input alm, with bouds (0:n-1,0:lmax,0:lmax)
!*    :itern (int) : Number of interation
!*
!* Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: fratio
  integer, intent(in) :: n, npix, lmax, chn
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: fratio = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(n,0:npix-1) :: nij
  double complex, intent(in), dimension(n,0:lmax,0:lmax) :: alm
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mgchain) :: mgc
  integer :: ni, l, ilmaxs(chn), insides(chn)
  double precision, dimension(n,0:lmax,0:lmax) :: clh
  double complex :: b(n,0:lmax,0:lmax)
  double precision :: ratio(itns(1))

  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides = (/int(dsqrt(npix/12d0))/)
  else
    ilmaxs  = lmaxs
    insides = nsides
  end if

  clh = 0d0
  do ni = 1, n
    do l = 1, lmax
      clh(ni,l,0:l) = dsqrt(cl(ni,l))
    end do
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi(n,npix,lmax,clh,nij,alm,b,'rhs')

  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call set_mgchain(mgc,chn,ilmaxs,insides,itns,eps,verbose)
  call cg_algorithm(n,clh,nij,b,xlm,mgc,1,ratio)
  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  !convert
  do l = 1, lmax
    if (filter=='')   xlm(:,l,0:l) = xlm(:,l,0:l)/clh(:,l,0:l)
    if (filter=='W')  xlm(:,l,0:l) = xlm(:,l,0:l)*clh(:,l,0:l)
  end do

end subroutine cnfiltertp


subroutine cnfilter_freq(n,fn,npix,lmax,cl,bl,nij,map,xlm,chn,lmaxs,nsides,itns,eps,fratio,filter,verbose)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is bxS+N
!*
!* Args:
!*    :n (int) : Number of maps
!*    :fn (int) : Maximum multipole of alm
!*    :npix (int) : Number of pixel
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :bl[fn,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :nij[n,fn,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:npix-1)
!*    :map[n,fn,pix] (double) : Input alm, with bouds (0:n-1,0:lmax,0:lmax)
!*
!* Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  logical, intent(in) :: verbose
  character(1), intent(in) :: filter
  character(100), intent(in) :: fratio
  integer, intent(in) :: n, fn, npix, lmax, chn
  integer, intent(in), dimension(1:chn) :: lmaxs, nsides, itns
  double precision, intent(in), dimension(1:chn) :: eps
  !opt4py :: chn = 1
  !opt4py :: lmaxs = [0]
  !opt4py :: nsides = [0]
  !opt4py :: itns = [1]
  !opt4py :: eps = [1e-6]
  !opt4py :: filter = ''
  !opt4py :: fratio = ''
  !opt4py :: verbose = False
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(fn,0:lmax) :: bl
  double precision, intent(in), dimension(n,fn,0:npix-1) :: nij, map
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  type(mgchain) :: mgc
  integer :: ni, fi, l, ilmaxs(chn), insides(chn)
  double precision, dimension(n,fn,0:lmax,0:lmax) :: clh
  double complex :: b(n,0:lmax,0:lmax)
  double precision :: ratio(itns(1))

  if (chn==1) then
    ilmaxs  = (/lmax/)
    insides = (/int(dsqrt(npix/12d0))/)
  else
    ilmaxs  = lmaxs
    insides = nsides
  end if

  clh = 0d0
  do ni = 1, n
    do l = 1, lmax
      do fi = 1, fn
        clh(ni,fi,l,0:l) = dsqrt(cl(ni,l))*bl(fi,l)
      end do
    end do
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi_freq_rhs(n,fn,npix,lmax,clh,nij,map,b)

  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  !call cg_algorithm_freq(n,k1,lmax,clh1,npix1,nij1,b1,xlm,itern,eps)
  call set_mgchain(mgc,chn,ilmaxs,insides,itns,eps,verbose)
  call cg_algorithm_freq(n,fn,clh,nij,b,xlm,mgc,1,ratio)
  if (fratio/='')  call savetxt(fratio,ratio,ow=.true.)
  call free_mgchain(mgc)

  do l = 1, lmax
    do ni = 1, n
      if (filter=='')   xlm(ni,l,0:l) = xlm(ni,l,0:l)/dsqrt(cl(ni,l))
      if (filter=='W')  xlm(ni,l,0:l) = xlm(ni,l,0:l)*dsqrt(cl(ni,l))
    end do
  end do

end subroutine cnfilter_freq


subroutine cnfilter_so(n,k1,k2,lmax,cl,bl1,bl2,npix1,npix2,nij1,nij2,map1,map2,itern,xlm,eps,filter)
!* Computing inverse-variance (default) or Wiener filtered multipoles: C^-1d
!* This code assumes
!*   1) The signal power spectrum is isotropic Gaussian. 
!*   2) Inverse noise covariance is given in pixel space and diagonal (nij = sigma x delta_ij).
!*   3) The data model is bxS+N
!*
!* Args:
!*    :n (int) : T(1), Q/U(2) or T/Q/U(3)
!*    :k1 (int) : Number of frequencies
!*    :k2 (int) : Number of frequencies
!*    :npix1 (int) : Number of pixels for each input maps and inv noise covariance
!*    :npix2 (int) : Number of pixels for each input maps and inv noise covariance
!*    :lmax (int) : Maximum multipole of alm
!*    :cl[n,l] (double) : Angular power spectrum of alm, with bounds (0:n-1,0:lmax)
!*    :bl1[k,l] (double) : Beam spectrum, with bounds (0:k1-1,0:lmax)
!*    :bl2[k,l] (double) : Beam spectrum, with bounds (0:k2-1,0:lmax)
!*    :nij1[n,k,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k1-1,0:npix1-1)
!*    :nij2[n,k,pix] (double) : Inverse of the noise variance at each pixel, with bounds (0:n-1,0:k2-1,0:npix2-1)
!*    :map1[n,k,pix] (double) : Input maps, with bouds (0:n-1,0:k1-1,0:npix1-1)
!*    :map2[n,k,pix] (double) : Input maps, with bouds (0:n-1,0:k2-1,0:npix2-1)
!*    :itern (int) : Number of interation
!*
!* Args(optional): 
!*    :eps (double): Numerical parameter to finish the iteration if ave(|Ax-b|)<eps, default to 1e-6
!*    :filter (str): C-inverse ('') or Wiener filter (W), default to C-inverse.
!*
!* Returns:
!*    :xlm[n,l,m] (dcmplx) : C-inverse / Wiener filtered multipoles, with bounds (0:n-1,0:lmax,0:lmax)
!*
  implicit none
  !I/O
  character(1), intent(in) :: filter
  integer, intent(in) :: n, k1, k2, npix1, npix2, lmax, itern
  double precision, intent(in) :: eps
  !opt4py :: eps = 1e-6
  !opt4py :: filter = ''
  double precision, intent(in), dimension(n,0:lmax) :: cl
  double precision, intent(in), dimension(k1,0:lmax) :: bl1
  double precision, intent(in), dimension(k2,0:lmax) :: bl2
  double precision, intent(in), dimension(n,k1,0:npix1-1) :: nij1, map1
  double precision, intent(in), dimension(n,k2,0:npix2-1) :: nij2, map2
  double complex, intent(out), dimension(n,0:lmax,0:lmax) :: xlm
  !internal
  integer :: ni, l, ki
  double precision, dimension(n,k1,0:lmax,0:lmax) :: clh1
  double precision, dimension(n,k2,0:lmax,0:lmax) :: clh2
  double complex, dimension(n,0:lmax,0:lmax) :: b1, b2

  clh1 = 0d0
  clh2 = 0d0
  do ni = 1, n
    do l = 1, lmax
      do ki = 1, k1
        clh1(ni,ki,l,0:l) = dsqrt(cl(ni,l))*bl1(ki,l)
      end do
      do ki = 1, k2
        clh2(ni,ki,l,0:l) = dsqrt(cl(ni,l))*bl2(ki,l)
      end do
    end do
  end do

  !first compute b = C^1/2 N^-1 X
  call mat_multi_freq_rhs(n,k1,npix1,lmax,clh1,nij1,map1,b1)
  call mat_multi_freq_rhs(n,k2,npix2,lmax,clh2,nij2,map2,b2)
  !solve x where [1 + C^1/2 N^-1 C^1/2] x = b
  call cg_algorithm_freq_nside(n,k1,k2,lmax,clh1,clh2,npix1,npix2,nij1,nij2,b1+b2,xlm,itern,eps)

  do l = 1, lmax
    do ni = 1, n
      if (filter=='')   xlm(ni,l,0:l) = xlm(ni,l,0:l)/dsqrt(cl(ni,l))
      if (filter=='W')  xlm(ni,l,0:l) = xlm(ni,l,0:l)*dsqrt(cl(ni,l))
    end do
  end do

end subroutine cnfilter_so


end module cninv
