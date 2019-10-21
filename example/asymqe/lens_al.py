# Flatsky lensing normalization from an asymmetric quadratic estimator
import numpy as np
import basic
import cmb

Tcmb  = 2.726e6    # CMB temperature
lmax  = 3000       # maximum multipole of output normalization
rlmin = 2
rlmax = 3000
#glmin = 500
#glmax = 2000
s0 = 70.
s1 = 20.
#t0 = .1
#t1 = .1

# load unlensed and lensed Cls
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2
TT  = lcl[0,:]
n0 = (s0*np.pi/10800./Tcmb)**2 #/ cmb.beam(t0,lmax)**2
n1 = (s1*np.pi/10800./Tcmb)**2 #/ cmb.beam(t1,lmax)**2
#n0[glmax+1:] = 1e30
#n1[:glmin+1] = 1e30

# compute analytic normalization
AA = TT + n0
BB = TT + n1
AB = TT
#np.savetxt('cl.dat',np.array((np.linspace(0,lmax,lmax+1),AA,BB)).T,fmt='%8.6e')

gln = 30
gle = 1e-10
print('PA')
Apa, Npa = basic.flat.alxy_asym('Tx','lensinga',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
Npa = Apa**2/Npa
print('AP')
Aap, Nap = basic.flat.alxy_asym('Tx','lensinga',lmax,rlmin,rlmax,TT,BB,AA,AB,gln,gle)
Nap = Aap**2/Nap
print('Comb')
Ac, Nc = basic.flat.alxy_asym('Tx','lensingc',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
#Ac, Nc = basic.flat.alxy_asym('Tx','lensingc',lmax,rlmin,rlmax,TT,BB,AA,AB,gln,gle)
Nx = Apa*Aap/Nc
Nc = (Npa*Nap-Nx**2)/(Npa+Nap-2*Nx)
'''
print('opt-0')
A0, N0 = basic.flat.alxy_asym('Tx','lensing0',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
N0 = A0**2/N0
print('opt-1')
A1, N1 = basic.flat.alxy_asym('Tx','lensing1',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
N1 = A1**2/N1
print('opt-2')
A2, N2 = basic.flat.alxy_asym('Tx','lensing2',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
N2 = A2**2/N2
print('opt-e')
Ae, Ne = basic.flat.alxy_asym('Tx','lensinge',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
Ne = Ae**2/Ne
print('opt')
Ao, No = basic.flat.alxy_asym('Tx','lensing',lmax,rlmin,rlmax,TT,AA,BB,AB,gln,gle)
'''

# save
#np.savetxt('al.dat',np.array((np.linspace(0,lmax,lmax+1),Npa,Nap,Nc,N0,N1,N2,Ne,Ao)).T,fmt='%8.6e')
np.savetxt('al.dat',np.array((np.linspace(0,lmax,lmax+1),Nx,Nc)).T,fmt='%8.6e')

