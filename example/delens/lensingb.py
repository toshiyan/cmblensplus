# * Linear template delensing

# load modules
import numpy as np
import basic
import curvedsky

# define parameters
Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output cl
rL   = [100,3000] # E-phi convolution multipole range

# load unlensed and lensed Cls
ucl  = basic.apstool.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2
lcl  = basic.apstool.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2

# generate gaussian CMB alms
Talm, Ealm, Balm = curvedsky.utils.sim_teb(lmax,lcl)

# generate gaussian phi
glm = curvedsky.utils.sim_alm(lmax,ucl[3,:])

# template lensing B-mode
lalm = curvedsky.delens.lensb(2048,rL,rL,lmax,Ealm,glm)

# aps
bb = basic.apstool.alm2cl(lmax,lalm,lalm)
gg = basic.apstool.alm2cl(lmax,glm,glm)

# save
np.savetxt('clbb.dat',np.array((np.linspace(0,lmax,lmax+1),bb,gg)).T)

