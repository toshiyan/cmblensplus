# A quick demonstration of lensing reconstruction in fullsky
# - this code shows how to compute the lensing estiamtors, and output their power spectra

# load modules
import numpy as np
import basic
import curvedsky

# define parameters
Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output normalization
rlmin = 100
rlmax = 3000      # reconstruction multipole range

# load unlensed and lensed Cls
#ucl  = basic.aps.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2
lcl  = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2

# calculate normalization (gradient and curl modes)
Ag = {}
Ac = {}
Ag['TT'], Ac['TT'] = curvedsky.norm_lens.qtt(lmax,rlmin,rlmax,lcl[0,:],lcl[0,:])
Ag['TE'], Ac['TE'] = curvedsky.norm_lens.qte(lmax,rlmin,rlmax,lcl[3,:],lcl[0,:],lcl[1,:])
Ag['TB'], Ac['TB'] = curvedsky.norm_lens.qtb(lmax,rlmin,rlmax,lcl[3,:],lcl[0,:],lcl[2,:])
Ag['EE'], Ac['EE'] = curvedsky.norm_lens.qee(lmax,rlmin,rlmax,lcl[1,:],lcl[1,:])
Ag['EB'], Ac['EB'] = curvedsky.norm_lens.qeb(lmax,rlmin,rlmax,lcl[1,:],lcl[1,:],lcl[2,:])

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
glm = {}
clm = {}
glm['TT'], clm['TT'] = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,lcl[0,:],Talm,Talm)
glm['TE'], clm['TE'] = curvedsky.rec_lens.qte(lmax,rlmin,rlmax,lcl[3,:],Talm,Ealm)
glm['TB'], clm['TB'] = curvedsky.rec_lens.qtb(lmax,rlmin,rlmax,lcl[3,:],Talm,Balm)
glm['EE'], clm['EE'] = curvedsky.rec_lens.qee(lmax,rlmin,rlmax,lcl[1,:],Ealm,Ealm)
glm['EB'], clm['EB'] = curvedsky.rec_lens.qeb(lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm)

# normalized estimators
for qest in ['TT','TE','TB','EE','EB']:
  glm[qest] *= Ag[qest][:,None]
  clm[qest] *= Ac[qest][:,None]

# compute power spectrum
gg = {}
cc = {}
for qest in ['TT','TE','TB','EE','EB']:
  gg[qest] = curvedsky.utils.alm2cl(lmax,glm[qest])
  cc[qest] = curvedsky.utils.alm2cl(lmax,clm[qest])

# save
np.savetxt('cl_grad.dat',np.array((np.linspace(0,lmax,lmax+1),gg['TT'],gg['TE'],gg['TB'],gg['EE'],gg['EB'])).T)
np.savetxt('cl_curl.dat',np.array((np.linspace(0,lmax,lmax+1),cc['TT'],cc['TE'],cc['TB'],cc['EE'],cc['EB'])).T)

