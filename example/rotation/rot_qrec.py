# A quick demonstration of cosmic birefringence reconstruction in fullsky
# - this code shows how to compute the estiamtors, and output their power spectra

# load modules
import numpy as np
import basic
import curvedsky

# define parameters
Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output normalization
rlmin = 100
rlmax = 3000      # reconstruction multipole range
L = np.linspace(0,lmax,lmax+1)

# load lensed Cls
lcl  = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2

# calculate normalization (gradient and curl modes)
Aa = {}
Aa['EB'] = curvedsky.norm_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],lcl[1,:],lcl[2,:])
np.savetxt('al.dat',np.array((L,Aa['EB'])).T,fmt='%4.6e')

# simple diagonal c-inverse
Fl = np.zeros((3,lmax+1,lmax+1))
for l in range(rlmin,rlmax):
  Fl[:,l,0:l+1] = 1./lcl[:3,l,None]

# load or generate CMB alms 
# - gaussian alms are generated, and the reconstructed cls are equal to the normalization
Talm, Ealm, Balm = curvedsky.utils.gaussTEB(lmax,lcl[0,:],lcl[1,:],lcl[2,:],lcl[3,:])
# - to measure lensing cl, you need Talm, Ealm and Balm from a lensed cmb simulation

# diagonal filtering (since idealistic)
Talm = Talm * Fl[0,:,:]
Ealm = Ealm * Fl[1,:,:]
Balm = Balm * Fl[2,:,:]

# compute unnormalized estiamtors
nsidet = 2048
alm = {}
alm['EB'] = curvedsky.rec_rot.qeb(lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm)

# normalized estimators
for qest in ['EB']:
  alm[qest] *= Aa[qest][:,None]

# compute power spectrum
aa = {}
for qest in ['EB']:
  aa[qest] = curvedsky.utils.alm2cl(lmax,alm[qest])

# save
np.savetxt('cl_rot.dat',np.array((L,aa['EB'])).T)

