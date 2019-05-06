# Flatsky lensing reconstruction normalization
# - This code compute normalization of estiamtors
import numpy as np
import flatsky
import basic

Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output normalization
rL   = [100,3000] # reconstruction multipole range
nx, ny  = 512, 512
D    = np.array([nx,ny])/60.*np.pi/180.
bn   = 50

oL   = [2,lmax]

# load unlensed and lensed Cls
lcl  = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
itt     = np.zeros(lmax+1)
iee     = np.zeros(lmax+1)
ibb     = np.zeros(lmax+1)
itt[2:] = 1./lcl[0,2:]
iee[2:] = 1./lcl[1,2:]
ibb[2:] = 1./lcl[2,2:]

# binned multipoles
bp, bc = basic.aps.binning(bn,oL)

# multipoles on grids
lx, ly, el, il = flatsky.utils.elarrays(nx,ny,D)
kl = el*(el+1.)/2.

# assign 1d cl on 2d grind
cltt = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,lcl[0,:])
clee = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,lcl[1,:])
clte = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,lcl[3,:])
iltt = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,itt)
ilee = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,iee)
ilbb = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,ibb)

# compute analytic normalization with 2d filtering
Ag = {}
Ac = {}
Ag['TT'], Ac['TT'] = flatsky.norm_lens.qtt(nx,ny,D,rL,iltt,cltt,oL)
Ag['TE'], Ac['TE'] = flatsky.norm_lens.qte(nx,ny,D,rL,iltt,ilee,clte,oL)
Ag['TB'], Ac['TB'] = flatsky.norm_lens.qtb(nx,ny,D,iltt,ilbb,clte,rL,oL)
Ag['EE'], Ac['EE'] = flatsky.norm_lens.qee(nx,ny,D,ilee,clee,rL,oL)
Ag['EB'], Ac['EB'] = flatsky.norm_lens.qeb(nx,ny,D,ilee,ilbb,clee,rL,oL)

# kappa binned spectrum
for q in ['TT','TE','TB','EE','EB']:
  Ag[q] = flatsky.utils.c2d2bcl(nx,ny,D,kl**2*Ag[q],bn,oL)
  Ac[q] = flatsky.utils.c2d2bcl(nx,ny,D,kl**2*Ac[q],bn,oL)

# save
np.savetxt('al_grad.dat',np.array((bc,Ag['TT'],Ag['TE'],Ag['TB'],Ag['EE'],Ag['EB'])).T,fmt='%8.6e')
np.savetxt('al_curl.dat',np.array((bc,Ac['TT'],Ac['TE'],Ac['TB'],Ac['EE'],Ac['EB'])).T,fmt='%8.6e')

