#!/usr/bin/env python
# coding: utf-8

import numpy as np
import cmblensplus.curvedsky as cs

# Generate random gaussian fields


# Set parameters

nside = 512  # Nside of Healpix map
npix  = 12*nside**2  # Numer of pixels of Healpix map
lmax  = 2*nside # Maximum multipole of harmonic coefficients to be computed

Tcmb = 2.72e6
nl = np.ones(lmax+1)*(10./Tcmb*np.pi/10800.)**2 # this is the white noise spectrum with 10 uK'

# Generate random gaussian fields
simn = 2  # nunmber of iterations
nmap = np.zeros((simn,npix)) # map for each random field
for i in range(simn):
    # generate random Gaussian fields from the angular power spectrum, nl
    nalm = cs.utils.gauss1alm(nl) # NOTE: alm is Healpix order, not healpy order
    # convert to Healpix map
    nmap[i,:] = cs.utils.hp_alm2map(nside,nalm)


cl = np.zeros((simn,lmax+1))
for i in range(simn):
    nalm = cs.utils.hp_map2alm(lmax,lmax,nmap[i,:])
    print('0,0',nalm[0,0],'1,0',nalm[1,0])
    cl[i,:] = cs.utils.alm2cl(nalm)
    print(np.min(cl[i,:]),np.max(cl[i,:]),np.shape(cl[i,:]))
    print(np.mean(cl[i,:]/nl))
np.savetxt('test_cl.dat',cl.T)


N = 3
Cov = np.zeros((N,N,lmax+1))
cl  = np.zeros((4,lmax+1))


# In[8]:


cl[:,2:lmax+1] = np.loadtxt('../data/lensedcls.dat',unpack=True,usecols=(1,2,3,4))[:,:lmax-1]
print(np.shape(cl))


# In[9]:


Cov[0,0,:] = cl[0,:lmax+1]
Cov[1,1,:] = cl[1,:lmax+1]
Cov[2,2,:] = cl[2,:lmax+1]
Cov[0,1,:] = Cov[1,0,:] = cl[3,:lmax+1]


# In[16]:


simn = 2  # nunmber of iterations
cls = np.zeros((4,simn,lmax+1)) # map for each random field
for i in range(simn):
    # generate random Gaussian fields from the angular power spectrum, nl
    alm = cs.utils.gaussalm(Cov)
    # compute cl
    for j in range(N):
        cls[j,i,:] = cs.utils.alm2cl(lmax,alm[j])
    cls[3,i,:] = cs.utils.alm2cl(lmax,alm[0],alm[1])
mcl = np.mean(cls,axis=1)


# In[22]:


plt.plot(cls[2,0,:])
plt.plot(Cov[1,0,:])


# In[ ]:




