# Fullsky cosmic birefringence reconstruction normalization
import numpy as np
import curvedsky
import basic

Tcmb  = 2.726e6    # CMB temperature
lmax  = 3000       # maximum multipole of output normalization
rlmin = 100        # reconstruction multipole range
rlmax = 3000

# load unlensed and lensed Cls
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
ocl = lcl

# calculate normalization (gradient and curl modes)
Ag = np.zeros((6,lmax+1))
Ag[4,:] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

# save
L = np.linspace(0,lmax,lmax+1)
np.savetxt('al.dat',np.concatenate((L[None,:],Ag)).T,fmt='%4.6e')

