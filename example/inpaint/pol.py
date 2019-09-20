# A quick demonstration of inpainting

import numpy as np
import healpy as hp
import basic
import curvedsky
from matplotlib.pyplot import *

# define parameters
Tcmb  = 2.726e6    # CMB temperature
lmax  = 4000       # maximum multipole of alm
nside = 4096
npix  = 12*nside**2
sigma = 3.         # uK'
L = np.linspace(0,lmax,lmax+1)

# load cl
cl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2

# load mask
mask = hp.fitsfunc.read_map('/project/projectdirs/sobs/delensing/mask/la.fits')
#mask = hp.fitsfunc.read_map('../data/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits')
mask = hp.pixelfunc.ud_grade(mask,nside)
w2 = np.average(mask**2)
print(w2)

# Gaussian alms are generated here
Talm, Ealm, Balm = curvedsky.utils.gaussTEB(lmax,cl[0,:],cl[1,:],cl[2,:],cl[3,:])

Tmap = mask * curvedsky.utils.hp_alm2map(npix,lmax,lmax,Talm)
Talm = curvedsky.utils.hp_map2alm(nside,lmax,lmax,Tmap)
Qmap, Umap = curvedsky.utils.hp_alm2map_spin(npix,lmax,lmax,2,Ealm,Balm)
Qmap *= mask
Umap *= mask
Ealm, Balm = curvedsky.utils.hp_map2alm_spin(nside,lmax,lmax,2,Qmap,Umap)

xl0  = curvedsky.utils.alm2cl(lmax,Talm)/w2
xl1  = curvedsky.utils.alm2cl(lmax,Ealm)/w2
xl2  = curvedsky.utils.alm2cl(lmax,Balm)/w2

# noise covariance
nij = mask * (sigma*(np.pi/10800.)/Tcmb)**(-2)
#nij = np.ones(npix)

# Wiener filtering
Talm, Ealm, Balm = curvedsky.cninv.cnfilterpol(3,npix,lmax,cl[0:3,:],np.array((nij,nij,nij)),np.array((Talm,Ealm,Balm)),2000,filter='W')

# power spectrum
Xl0 = curvedsky.utils.alm2cl(lmax,Talm)/w2
Xl1 = curvedsky.utils.alm2cl(lmax,Ealm)/w2
Xl2 = curvedsky.utils.alm2cl(lmax,Balm)/w2

# filtered map
#xmap = curvedsky.utils.hp_alm2map(npix,lmax,lmax,xlm)
#hp.mollview(xmap)
#hp.mollview(xmap-Tmap)
#show()
np.savetxt('cl.dat',np.array((L,cl[0,:],xl0,Xl0,xl1,Xl1,xl2,Xl2)).T)

