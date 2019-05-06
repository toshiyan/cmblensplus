# Fullsky lensing reconstruction normalization
# - This code compute normalization of estiamtors
import numpy as np
import curvedsky
import basic

Tcmb  = 2.726e6    # CMB temperature
lmax  = 1000       # maximum multipole of output normalization
rlmin = 100        # reconstruction multipole range
rlmax = 1000
QDO = [True,True,True,True,True,False]
QDO = [True,True,True,False,False,False]

# load unlensed and lensed Cls
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
ocl = lcl

'''
# calculate normalization (gradient and curl modes)
Ag = np.zeros((6,lmax+1))
Ac = np.zeros((6,lmax+1))
Ig = np.zeros((4,lmax+1))
Ic = np.zeros((4,lmax+1))
Ag[0,:], Ac[0,:] = curvedsky.norm_lens.qtt(lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
Ag[1,:], Ac[1,:] = curvedsky.norm_lens.qte(lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[1,:])
Ag[2,:], Ac[2,:] = curvedsky.norm_lens.qee(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:])
Ag[3,:], Ac[3,:] = curvedsky.norm_lens.qtb(lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:])
Ag[4,:], Ac[4,:] = curvedsky.norm_lens.qeb(lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

# mv
Ig[0,:], Ic[0,:] = curvedsky.norm_lens.qttte(lmax,rlmin,rlmax,lcl[0,:],lcl[3,:],ocl[0,:],ocl[1,:],ocl[3,:])
Ig[1,:], Ic[1,:] = curvedsky.norm_lens.qttee(lmax,rlmin,rlmax,lcl[0,:],lcl[1,:],ocl[0,:],ocl[1,:],ocl[3,:])
Ig[2,:], Ic[2,:] = curvedsky.norm_lens.qteee(lmax,rlmin,rlmax,lcl[1,:],lcl[3,:],ocl[0,:],ocl[1,:],ocl[3,:])
Ig[3,:], Ic[3,:] = curvedsky.norm_lens.qtbeb(lmax,rlmin,rlmax,lcl[1,:],lcl[2,:],lcl[3,:],ocl[0,:],ocl[1,:],ocl[2,:],ocl[3,:])
Ag[5,:], nlg = curvedsky.norm_lens.qmv(lmax,QDO,Ag,Ig)
Ac[5,:], nlc = curvedsky.norm_lens.qmv(lmax,QDO,Ac,Ic)
'''

Ag, Ac, nlg, nlc = curvedsky.norm_lens.qall(QDO,lmax,rlmin,rlmax,lcl,ocl)

# save
L = np.linspace(0,lmax,lmax+1)
np.savetxt('al.dat',np.concatenate((L[None,:],Ag,Ac)).T,fmt='%4.6e')
#np.savetxt('il.dat',np.concatenate((L[None,:],Ig,Ic)).T,fmt='%4.6e')
np.savetxt('nl.dat',np.concatenate((L[None,:],nlg,nlc)).T,fmt='%4.6e')

