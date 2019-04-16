# A simple lensing reconstruction in flat sky
import numpy as np
import flatskypy as flatsky
import basicpy as basic

# parameters
Tcmb = 2.72e6
lmax = 3000
rL = [500,3000]
oL = [2,3000]
nx, ny = 512, 512
D  = np.array([nx,ny]) / 120.*np.pi/180.
bn = 50
mc = 10

# binned multipoles
bp, bc = basic.aps.binning(bn,oL)

# multipoles on grids
lx, ly, el, il = flatsky.utils.elarrays(nx,ny,D)
kl = el*(el+1.)/2.

# load unlensed and lensed Cls
lcl  = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
itt     = np.zeros(lmax+1)
iee     = np.zeros(lmax+1)
ibb     = np.zeros(lmax+1)
itt[2:] = 1./lcl[0,2:]
iee[2:] = 1./lcl[1,2:]
ibb[2:] = 1./lcl[2,2:]

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

# loop for realizations
cls = np.zeros((mc,3,bn))
cks = np.zeros((mc,5,bn))

# loop over MC realizations
for i in range(mc):

  print(i)
  # generate Fourier mode on 2d grids
  tlm, elm = flatsky.utils.gauss2alm(nx,ny,D,2,lmax,lcl[0,:],lcl[3,:],lcl[1,:])
  blm = flatsky.utils.gauss1alm(nx,ny,D,2,lmax,lcl[2,:])

  #cls[i,0,:] = flatsky.utils.alm2bcl(bn,oL,nx,ny,D,tlm)
  #cls[i,1,:] = flatsky.utils.alm2bcl(bn,oL,nx,ny,D,tlm,elm)
  #cls[i,2,:] = flatsky.utils.alm2bcl(bn,oL,nx,ny,D,elm)

  # filtering
  tlm *= iltt
  elm *= ilee
  blm *= ilbb

  # reconstruction
  klm = {}
  klm['TT'], clm = flatsky.rec_lens.qtt(nx,ny,D,rL,cltt,tlm,tlm,gtype='k')
  klm['TE'], clm = flatsky.rec_lens.qte(nx,ny,D,rL,clte,tlm,elm,gtype='k')
  klm['TB'], clm = flatsky.rec_lens.qtb(nx,ny,D,rL,clte,tlm,blm,gtype='k')
  klm['EE'], clm = flatsky.rec_lens.qee(nx,ny,D,rL,clee,elm,elm,gtype='k')
  klm['EB'], clm = flatsky.rec_lens.qeb(nx,ny,D,rL,clee,elm,blm,gtype='k')

  for I, q in enumerate(['TT','TE','TB','EE','EB']):
    klm[q] *= Ag[q]*kl**2
    cks[i,I,:] = flatsky.utils.alm2bcl(bn,oL,nx,ny,D,klm[q])

# save
#mcls = np.mean(cls,axis=0)
mcks = np.mean(cks,axis=0)
#np.savetxt('cl_cmb.dat',np.array((bc,mcls[0,:],mcls[1,:],mcls[2,:])).T)
np.savetxt('cl_krec.dat',np.concatenate((bc[None,:],mcks)).T)

