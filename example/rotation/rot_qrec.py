# A quick demonstration of cosmic birefringence reconstruction in fullsky
# - this code shows how to compute the estiamtors, and output their power spectra

# load modules
import numpy as np
import basic
import curvedsky

# define parameters
Tcmb  = 2.726e6    # CMB temperature
lmax  = 2048       # maximum multipole of output normalization
rlmin = 100
rlmax = 2048      # reconstruction multipole range
sig   = 10.
adeg  = 1. # signal
ac2rad = np.pi/180./60.
L = np.linspace(0,lmax,lmax+1)

# load unlensed and lensed Cls
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
nl  = np.zeros((4,lmax+1))
nl[0,:] = (sig*ac2rad/Tcmb)**2
nl[1,:] = 2*nl[0,:]
nl[2,:] = 2*nl[0,:]
ocl = lcl + nl

# calculate normalizations
Al = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

# simple diagonal c-inverse
Fl = np.zeros((3,lmax+1,lmax+1))
for l in range(rlmin,rlmax):
    Fl[:,l,0:l+1] = 1./ocl[:3,l,None]

# generate CMB alms 
# - gaussian alms are generated, and the reconstructed cls are equal to the normalization
Talm, Ealm, Balm = curvedsky.utils.gaussTEB(lmax,lcl[0,:],lcl[1,:],lcl[2,:],lcl[3,:])
Talm += curvedsky.utils.gauss1alm(lmax,nl[0,:])
Ealm += curvedsky.utils.gauss1alm(lmax,nl[1,:])
Balm += curvedsky.utils.gauss1alm(lmax,nl[2,:])
aalm = curvedsky.utils.gauss1alm(lmax,np.exp(-L**2*(adeg*np.pi/180.)**2))

alpha = curvedsky.utils.hp_alm2map(12*2048**2,lmax,lmax,aalm)
Q, U = curvedsky.utils.hp_alm2map_spin(12*2048**2,lmax,lmax,2,Ealm,Balm)
rQ = Q*np.cos(2.*alpha) - U*np.sin(2.*alpha)
rU = Q*np.sin(2.*alpha) + U*np.cos(2.*alpha)
Ealm, Balm = curvedsky.utils.hp_map2alm_spin(2048,lmax,lmax,2,rQ,rU)

# diagonal filtering
Talm *= Fl[0,:,:]
Ealm *= Fl[1,:,:]
Balm *= Fl[2,:,:]

# compute unnormalized estiamtors
alm = curvedsky.rec_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm)

# normalized estimators
alm *= Al[:,None]

# compute power spectrum
cl = curvedsky.utils.alm2cl(lmax,alm,aalm)

# save
np.savetxt('cl_rot.dat',np.array((L,cl)).T)

