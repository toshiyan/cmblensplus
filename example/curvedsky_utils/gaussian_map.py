#!/usr/bin/env python
# coding: utf-8

# In[1]:


# external
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
# modules from curvedsky/wrap/
import cmblensplus.basic as basic
import cmblensplus.curvedsky as cs
import constant as c
import cmb


# Set parameters

# In[2]:


nside = 1024  # Nside of Healpix map
npix  = 12*nside**2  # Numer of pixels of Healpix map
lmax  = 2*nside # Maximum multipole of harmonic coefficients to be computed


# In[3]:


# ucl is an array of shape [0:5,0:rlmax+1] and ucl[0,:] = TT, ucl[1,:] = EE, ucl[2,:] = TE, lcl[3,:] = phiphi, lcl[4,:] = Tphi
ucl = cmb.read_camb_cls('../data/unlensedcls.dat',ftype='scal',output='array')[:,:lmax+1]


# Generate random gaussian fields

# In[4]:


Talm = cs.utils.gauss1alm(lmax,ucl[0,:])
Tmap = cs.utils.hp_alm2map(nside,lmax,lmax,Talm)


# In[5]:


palm = cs.utils.gauss1alm(lmax,ucl[3,:])
pmap = cs.utils.hp_alm2map(nside,lmax,lmax,palm)


# In[6]:


hp.mollview(Tmap)


# In[7]:


hp.mollview(Tmap*pmap)


# In[ ]:




