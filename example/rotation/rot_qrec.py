# A quick demonstration of cosmic birefringence reconstruction in fullsky
# - this code shows how to compute the estiamtors, and output their power spectra

# load modules
import numpy as np
import basic
import curvedsky

# define parameters
Tcmb  = 2.726e6    # CMB temperature
lmax  = 3000       # maximum multipole of output normalization
rlmin = 100
rlmax = 3000      # reconstruction multipole range
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

# calculate normalizations
Al = {}
Rl = {}
Al['qEB'] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])
Al['oEB'] = curvedsky.norm_tau.oeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])
Rl['qEB'] = curvedsky.norm_rot.teb(lmax,rlmin,rlmax,lcl[1,:],lcl[1,:],ocl[1,:],ocl[2,:])
Il = 1./(1.-Al['qEB']*Al['oEB']*Rl['qEB']**2)
Al['QEB'] = Al['qEB']*Il
np.savetxt('al.dat',np.array((L,Al['qEB'],Al['oEB'],Rl['qEB'],Al['QEB'])).T,fmt='%4.6e')

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

# diagonal filtering
Talm *= Fl[0,:,:]
Ealm *= Fl[1,:,:]
Balm *= Fl[2,:,:]

# compute unnormalized estiamtors
nsidet = 2048
alm = {}
alm['qEB'] = curvedsky.rec_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm)
alm['oEB'] = curvedsky.rec_tau.oeb(lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm)

# normalized estimators
for qest in ['qEB','oEB']:
    alm[qest] *= Al[qest][:,None]

# bhe
alm['QEB'] = alm['qEB'] - Al['qEB'][:,None]*Rl['qEB'][:,None]*alm['oEB']
alm['QEB'] *= Il[:,None]

# compute power spectrum
cl = {}
for qest in ['qEB','oEB','QEB']:
    cl[qest] = curvedsky.utils.alm2cl(lmax,alm[qest])

# save
np.savetxt('cl_rot.dat',np.array((L,cl['qEB'],cl['oEB'],cl['QEB'])).T)

