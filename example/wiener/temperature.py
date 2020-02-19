#!/usr/bin/env python
# coding: utf-8

# ### A quick demonstration of temperature wiener filtering

import numpy as np
import healpy as hp
from matplotlib.pyplot import *
# from cmblensplus
import basic
import curvedsky


# define parameters
Tcmb  = 2.726e6 # CMB temperature
lmax  = 1000     # maximum multipole of alm
#lmax = 8
#lmax  = 400
nside = 512
#nside = 4
#nside = 256
npix  = 12*nside**2
sigma = 100.    # uK'
L = np.linspace(0,lmax,lmax+1)


# load input cl
cl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2


# load mask
mask = hp.fitsfunc.read_map('../data/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits',verbose=False)
mask = hp.pixelfunc.ud_grade(mask,nside)
w2 = np.average(mask**2)
print('W2 correction is',w2)


# generate Temperature fluctuations
Talm = curvedsky.utils.gauss1alm(lmax,cl[0,:])


# CMB temperature map with mask
Tmap = mask * curvedsky.utils.hp_alm2map(npix,lmax,lmax,Talm)
Talm = curvedsky.utils.hp_map2alm(nside,lmax,lmax,Tmap)

#cl0 = curvedsky.utils.alm2cl(lmax,Talm)/w2
#print(cl[0,2:50]/(sigma*(np.pi/10800.)/Tcmb)**2)

# inverse noise covariance
#nij  = np.ones(npix)*(sigma*(np.pi/10800.)/Tcmb)**(-2)
nij  = mask * (sigma*(np.pi/10800.)/Tcmb)**(-2)


# Wiener filtering
import time
t1 = time.time()
xlm0 = curvedsky.cninv.cnfilter(npix,lmax,cl[0,:],nij,Talm,itns=[300],eps=[1e-6],fratio='cnv_fid.dat',filter='W')
t2 = time.time()
print(t2-t1)

# Wiener filtering with multigrid preconditioning
'''
eps = [1e-6,0.]
itns = [100,0]
#lmaxs = [lmax,16]
lmaxs = [lmax,32]
#nsides = [nside,nside]
nsides = [nside,64]
t1 = time.time()
xlm1 = curvedsky.cninv.cnfilter(npix,lmax,cl[0,:],nij,Talm,2,lmaxs,nsides,itns,eps,filter='W')
t2 = time.time()
print(t2-t1)
'''

#'''
eps = [1e-6,.1,.1,.1,0.]
for itnmax in [1,2,3,4,5,6,7]:
    for lsp in [20,30,40,50,60]:
        t1 = time.time()
        lmaxs = [lmax,400,200,100,lsp]
        nsides = [nside,256,128,64,64]
        itns = [100,7,5,itnmax,0]
        f = 'cnv_n4_l'+str(lsp)+'_i'+str(itnmax)+'.dat'
        xlm1 = curvedsky.cninv.cnfilter(npix,lmax,cl[0,:],nij,Talm,5,lmaxs,nsides,itns,eps,filter='W',verbose=False,fratio=f)
        t2 = time.time()
        print(t2-t1)
#'''

'''
eps = [1e-6,.1,.1,.1,0.]
itns = [100,10,10,5,0]
lmaxs = [lmax,400,200,100,8]
nsides = [nside,256,128,64,16]
xlm1 = curvedsky.cninv.cnfilter(npix,lmax,cl[0,:],nij,Talm,5,lmaxs,nsides,itns,eps,filter='W')
'''

'''
eps = [1e-6,0.1,0.1,0.]
#for lmax0 in [600,500,400,300]:
for lmax0 in [200,150,100]:
    for itnmax in [1,2,3,4,5,6,7,8,9,10]:
        itns = [100,8,itnmax,0]
        lmaxs = [lmax,500,lmax0,8]
        nsides = [nside,256,128,128]
        t1 = time.time()
        f = 'cnv_n3_l'+str(lmax0)+'_i'+str(itnmax)+'.dat'
        xlm1 = curvedsky.cninv.cnfilter(npix,lmax,cl[0,:],nij,Talm,4,lmaxs,nsides,itns,eps,filter='W',fratio=f)
        t2 = time.time()
        print(t2-t1)
'''

