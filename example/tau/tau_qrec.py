# A quick demonstration of amplitude modulation reconstruction in fullsky
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
Al['tTT'] = curvedsky.norm_tau.qtt(lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
Al['lTT'], c = curvedsky.norm_lens.qtt(lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
#Al['tEB'] = curvedsky.norm_tau.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])
#Al['lEB'], c = curvedsky.norm_lens.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])
Al['xTT'] = curvedsky.norm_lens.ttt(lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
Il = 1./(1.-Al['tTT']*Al['lTT']*Al['xTT']**2)
Al['qTT'] = Al['tTT']*Il
#Il = 1./(1.-Al['lEB']*Al['tEB']*Al['xEB']**2)
#Al['lEB'] = Al['lEB']*Il
np.savetxt('al.dat',np.array((L,Al['tTT'],Al['qTT'])).T,fmt='%4.6e')

# simple diagonal c-inverse
Fl = np.zeros((3,lmax+1,lmax+1))
for l in range(rlmin,rlmax):
    Fl[:,l,0:l+1] = 1./ocl[:3,l,None]

# generate CMB alms 
# - gaussian alms are generated, and the reconstructed cls are equal to the normalization
Talm, Ealm, Balm = curvedsky.utils.gaussTEB(lmax,lcl[0,:],lcl[1,:],lcl[2,:],lcl[3,:])
tlm = curvedsky.utils.gauss1alm(lmax,1e-5*np.exp(-(L/500.)**2))

# modulate amplitude in map space
nside = 2048
npix  = 12*nside**2
tau  = curvedsky.utils.hp_alm2map(npix,lmax,lmax,tlm)
Tmap = curvedsky.utils.hp_alm2map(npix,lmax,lmax,Talm)
Talm = curvedsky.utils.hp_map2alm(nside,lmax,lmax,Tmap*(1.+tau))
Qmap, Umap = curvedsky.utils.hp_alm2map_spin(npix,lmax,lmax,2,Ealm,Balm)
Ealm, Balm = curvedsky.utils.hp_map2alm_spin(nside,lmax,lmax,2,Qmap*(1.+tau),Umap*(1.+tau))

# add noise
Talm += curvedsky.utils.gauss1alm(lmax,nl[0,:])
Ealm += curvedsky.utils.gauss1alm(lmax,nl[1,:])
Balm += curvedsky.utils.gauss1alm(lmax,nl[2,:])

# diagonal filtering
Talm *= Fl[0,:,:]
Ealm *= Fl[1,:,:]
Balm *= Fl[2,:,:]

# compute unnormalized estiamtors
alm = {}
alm['tTT'] = curvedsky.rec_tau.qtt(lmax,rlmin,rlmax,lcl[0,:],Talm,Talm)
#alm['qEB'] = curvedsky.rec_tau.qeb(lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm)
alm['lTT'], clm = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,lcl[0,:],Talm,Talm)

# normalized estimators
for qest in ['lTT','tTT']:
    alm[qest] *= Al[qest][:,None]

# bhe
alm['qTT'] = alm['tTT'] - Al['tTT'][:,None]*Al['xTT'][:,None]*alm['lTT']
alm['qTT'] *= Il[:,None]

# compute auto/cross spectra
cl = {}
xl = {}
for qest in ['qTT','tTT']:
    cl[qest] = curvedsky.utils.alm2cl(lmax,alm[qest])
    xl[qest] = curvedsky.utils.alm2cl(lmax,alm[qest],tlm)
tl = curvedsky.utils.alm2cl(lmax,tlm)
al = curvedsky.utils.alm2cl(lmax,alm['lTT'],tlm)

# save
np.savetxt('cl.dat',np.array((L,cl['tTT'],cl['qTT'],xl['tTT'],xl['qTT'],tl,al)).T)

