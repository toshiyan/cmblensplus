# Fullsky cosmic birefringence reconstruction normalization
import numpy as np
import curvedsky
import basic

Tcmb  = 2.726e6    # CMB temperature
lmax  = 3000       # maximum multipole of output normalization
rlmin = 300        # reconstruction multipole range
rlmax = 3000
sig   = 10.
ac2rad = np.pi/180./60.
L = np.linspace(0,lmax,lmax+1)

# load unlensed and lensed Cls
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
nl  = np.zeros((4,lmax+1))
nl[0,:] = (sig*ac2rad/Tcmb)**2
nl[1,:] = 2*nl[0,:]
nl[2,:] = 2*nl[0,:]
ocl = lcl + nl
#eb = lcl[1,:]
#eb = (L/80.)**(-2.42)
eb = (L**(-5)/80.)

# calculate normalizations
Al = np.zeros((4,lmax+1))
#Al[0,:] = curvedsky.norm_rot.qtb(lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:])
Al[0,:] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])
#Al[1,:] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:],lcl[2,:])
Al[1,:] = curvedsky.norm_tau.oeb(lmax,rlmin,rlmax,eb,ocl[1,:],ocl[2,:])
Al[2,:] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,eb,ocl[1,:],ocl[2,:],-eb) #mask MF
Rl = curvedsky.norm_rot.teb(lmax,rlmin,rlmax,lcl[1,:],eb,ocl[1,:],ocl[2,:]) #alpha x mask
Al[3,:] = Al[0,:]/(1.-Rl**2*Al[0,:]*Al[1,:])

# save
np.savetxt('al.dat',np.concatenate((L[None,:],Al)).T,fmt='%4.6e')

