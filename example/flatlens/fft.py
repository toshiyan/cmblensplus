# DFT
# - This code compute normalization of estiamtors
import numpy as np
import flatsky_py2 as flatsky
import basictools_py2 as basictools

Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output normalization
nx   = 2580
ny   = 1501
D    = np.array([nx,ny])/30.*np.pi/180.
bn   = 10
oL   = [2,lmax]

# load unlensed and lensed Cls
ucl  = basictools.apstool.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2

# binned multipoles
bp, bc = basictools.apstool.binning(bn,oL)

# multipoles on grids
lx, ly, el, il = flatsky.utils.el2dgrids(nx,ny,D)

# assign 1d cl on 2d grind
cltt = basictools.apstool.cl_to_c2d(nx,ny,el,ucl[0,:],oL)

# compute analytic normalization with 2d filtering
dltt = flatsky.fft.dft2dr(cltt,nx,ny,D,-1)
iltt = flatsky.fft.dft2dr(dltt,nx,ny,D,1)

# binned spectrum
cl0 = basictools.apstool.c2d_to_bcl(bn,oL,el,cltt)
cl1 = basictools.apstool.c2d_to_bcl(bn,oL,el,iltt)

# save
np.savetxt('fft.dat',np.array((bc,cl0,cl1)).T)

