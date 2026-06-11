# Flatsky reconstruction normalization for tau and point src
import numpy as np
import flatsky
import basic

Tcmb = 2.726e6     # CMB temperature
lmax = 3000        # maximum multipole of output normalization
rL   = [100,3000] # reconstruction multipole range
nx, ny  = 1024, 1024
D    = np.array([nx,ny])/60.*np.pi/180.
bn   = 30

oL   = [2,lmax]

# load unlensed and lensed Cls
#lcl     = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
lcl     = basic.aps.read_cambcls('test.dat',2,lmax,4,bb=True)/Tcmb**2
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

# assign 1d cl on 2d grind
cltt = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,lcl[0,:])
clee = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,lcl[1,:])
clte = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,lcl[3,:])
iltt = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,itt)
ilee = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,iee)
ilbb = flatsky.utils.cl2c2d(nx,ny,D,2,lmax,ibb)

# compute analytic normalization with 2d filtering
cltt += 1e-15
A = {}
for q in ['TT']:
  A[q] = {}
A['TT']['kk'], A['TT']['cc'] = flatsky.norm_lens.qtt(nx,ny,D,rL,iltt,cltt,oL)
A['TT']['kk'] *= (el*(el+1.)/2.)**2
A['TT']['cc'] *= (el*(el+1.)/2.)**2
A['TT']['tt'] = flatsky.norm_tau.qtt(nx,ny,D,rL,iltt,cltt,oL)
A['TT']['ss'] = flatsky.norm_src.qtt(nx,ny,D,rL,iltt,oL)
A['TT']['kt'] = flatsky.norm_kxt.qtt(nx,ny,D,rL,iltt,cltt,oL)
A['TT']['ks'] = flatsky.norm_kxs.qtt(nx,ny,D,rL,iltt,cltt,oL)
#Ag['TE'], Ac['TE'] = flatsky.norm_lens.qte(nx,ny,D,rL,iltt,ilee,clte,oL)
#Ag['TB'], Ac['TB'] = flatsky.norm_lens.qtb(nx,ny,D,iltt,ilbb,clte,rL,oL)
#Ag['EE'], Ac['EE'] = flatsky.norm_lens.qee(nx,ny,D,ilee,clee,rL,oL)
#Ag['EB'], Ac['EB'] = flatsky.norm_lens.qeb(nx,ny,D,ilee,ilbb,clee,rL,oL)

# kappa binned spectrum
for q in ['TT']:
#for q in ['TT','TE','TB','EE','EB']:
  for mc in ['kk','cc','tt','ss','kt','ks']:
    A[q][mc] = flatsky.utils.c2d2bcl(nx,ny,D,A[q][mc],bn,oL)

    # save
    np.savetxt('al_'+mc+'.dat',np.array((bc,A[q][mc])).T,fmt='%8.6e')

