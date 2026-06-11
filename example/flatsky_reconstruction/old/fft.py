# DFT
# - This code compute normalization of estiamtors
import numpy as np
import flatsky
import basic
import cmb
import binning

Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output normalization
nx   = 258
ny   = 150
D    = np.array([nx,ny])/30.*np.pi/180.
bn   = 10
oL   = [2,lmax]

# load unlensed and lensed Cls
ucl = cmb.read_camb_cls('../data/unlensedcls.dat',ftype='scal',output='array')[:,:lmax+1]

# binned multipoles
bp, bc = binning.binned_ells(bn,oL[0],oL[1])

# multipoles on grids
lx, ly, el, il = flatsky.utils.elarrays(nx,ny,D)

# assign 1d cl on 2d grind
cltt = flatsky.utils.cl2c2d(nx,ny,D,oL[0],oL[1],ucl[0,:])

print('A')

dltt = flatsky.ffttools.dft2dr(cltt,nx,ny,D,-1)
print('A')
iltt = flatsky.ffttools.dft2dr(dltt,nx,ny,D,1)
print('A')

# binned spectrum
print('A')

cl0 = flatsky.utils.c2d2bcl(nx,ny,D,cltt,bn,oL)
print('A')

cl1 = flatsky.utils.c2d2bcl(nx,ny,D,iltt,bn,oL)

# save
np.savetxt('fft.dat',np.array((bc,cl0,cl1)).T)

