# DFT
# - This code compute normalization of estiamtors
import numpy as np
import flatsky
import basic

Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output normalization
nx   = 2580
ny   = 1501
D    = np.array([nx,ny])/30.*np.pi/180.
bn   = 10
oL   = [2,lmax]

# load unlensed and lensed Cls
ucl  = basic.aps.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2

# binned multipoles
bp, bc = basic.aps.binning(bn,oL)

# multipoles on grids
lx, ly, el, il = flatsky.utils.el2dgrids(nx,ny,D)

# assign 1d cl on 2d grind
cltt = flatsky.utils.cl2c2d(nx,ny,el,ucl[0,:],oL)

dltt = flatsky.fft.dft2dr(cltt,nx,ny,D,-1)
iltt = flatsky.fft.dft2dr(dltt,nx,ny,D,1)

# binned spectrum
cl0 = flatsky.utils.c2d2bcl(bn,oL,el,cltt)
cl1 = flatsky.utils.c2d2bcl(bn,oL,el,iltt)

# save
np.savetxt('fft.dat',np.array((bc,cl0,cl1)).T)

