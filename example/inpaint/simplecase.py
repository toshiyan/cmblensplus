# A quick demonstration of inpainting

import numpy as np
import curvedsky

# define parameters
lmax = 1000       # maximum multipole of alm
npix = 12*512**2

# Gaussian alms are generated here
cl  = np.ones(lmax+1)
clh = np.ones((lmax+1,lmax+1))
alm = curvedsky.utils.gauss1alm(lmax,cl)
nij = np.ones(npix)
xlm = curvedsky.cninv.cg_algorithm(npix,lmax,clh,nij,alm,1000)
xl = curvedsky.utils.alm2cl(lmax,xlm)
np.savetxt('test.dat',xl)

