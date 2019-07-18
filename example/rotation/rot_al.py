# Fullsky cosmic birefringence reconstruction normalization
import numpy as np
import curvedsky
import basic

Tcmb  = 2.726e6    # CMB temperature
lmax  = 3000       # maximum multipole of output normalization
rlmin = 400        # reconstruction multipole range
rlmax = 3000
sig   = 10.
ac2rad = np.pi/180./60.

# load unlensed and lensed Cls
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
nl  = np.zeros((4,lmax+1))
nl[0,:] = (sig*ac2rad/Tcmb)**2
nl[1,:] = 2*nl[0,:]
nl[2,:] = 2*nl[0,:]
ocl = lcl + nl
print(lcl[2,:]/nl[2,:])

# calculate normalization (gradient and curl modes)
Ag = np.zeros((6,lmax+1))
Ag[3,:] = curvedsky.norm_rot.qtb(lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:])
Ag[4,:] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

# save
L = np.linspace(0,lmax,lmax+1)
np.savetxt('al.dat',np.concatenate((L[None,:],Ag)).T,fmt='%4.6e')

