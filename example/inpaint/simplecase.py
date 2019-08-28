# A quick demonstration of inpainting

import numpy as np
import curvedsky

# define parameters
lmax = 1000       # maximum multipole of alm
npix = 12*512**2

# Gaussian alms are generated here
cl  = np.ones(lmax+1)
alm = curvedsky.utils.gauss1alm(lmax,cl)
nij = np.ones(npix)
xlm = curvedsky.utils.cg_algorithm(npix,lmax,cl,nij,alm,10)
xl = curvedsky.utils.alm2cl(lmax,xlm)
np.savetxt('test.dat',xl)

