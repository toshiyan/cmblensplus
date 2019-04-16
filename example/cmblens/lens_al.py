# Fullsky lensing reconstruction normalization
# - This code compute normalization of estiamtors
import numpy as np
import devanalyticpy as analytic
import basictools

Tcmb  = 2.726e6    # CMB temperature
lmax  = 3000       # maximum multipole of output normalization
rlmin = 100        # reconstruction multipole range
rlmax = 3000

# load unlensed and lensed Cls
ucl  = basictools.apstool.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2
lcl  = basictools.apstool.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2

# calculate normalization (gradient and curl modes)
gltt, cltt = analytic.rec_lens.qtt(lmax,rlmin,rlmax,lcl[0,:],lcl[0,:])
glte, clte = analytic.rec_lens.qte(lmax,rlmin,rlmax,lcl[3,:],lcl[0,:],lcl[1,:])
gltb, cltb = analytic.rec_lens.qtb(lmax,rlmin,rlmax,lcl[3,:],lcl[0,:],lcl[2,:])
glee, clee = analytic.rec_lens.qee(lmax,rlmin,rlmax,lcl[1,:],lcl[1,:])
gleb, cleb = analytic.rec_lens.qeb(lmax,rlmin,rlmax,lcl[1,:],lcl[1,:],lcl[2,:])

# save
np.savetxt('al_grad.dat',np.array((np.linspace(0,lmax,lmax+1),gltt,glte,gltb,glee,gleb)).T)
np.savetxt('al_curl.dat',np.array((np.linspace(0,lmax,lmax+1),cltt,clte,cltb,clee,cleb)).T)

