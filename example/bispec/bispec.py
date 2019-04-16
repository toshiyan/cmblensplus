# Example of computing binned reduced bispectrum

# load modules
import numpy as np
import basictools
import curvedsky

# define parameters
Tcmb = 2.726e6    # CMB temperature
lmax = 1000       # maximum multipole of output normalization
bn   = 10         # number of multipole bins

# choose multipoles within a multipole bin
print('define multipole bins')
bmin = np.array([np.int(lmax*(i/bn)) for i in range(bn)])
bmax = bmin + np.int(lmax/bn)
bc   = (bmin + bmax)*.5
bp   = np.array([np.int(lmax*(i/bn)) for i in range(bn+1)])
sL   = bp[:2]

# load unlensed and lensed Cls
ucl  = basictools.apstool.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2

# generate gaussian phi
glm = curvedsky.utils.sim_alm(lmax,ucl[3,:])

# compute pmap = gmap + gmap**2
print('quad gauss fields')
plm = curvedsky.bispec.quad_gauss(lmax,glm)

# compute binned bispectra
bl = np.zeros((4,bn))

print('equilateral')
hl = curvedsky.bispec.equi_binnorm(bn,bp)
bl[0,:] = curvedsky.bispec.equi_bin(bn,bp,plm) * np.sqrt(4*np.pi)/hl

print('fold')
hl      = curvedsky.bispec.fold_binnorm(bn,bp)
bl[1,:] = curvedsky.bispec.fold_bin(bn,bp,plm) * np.sqrt(4*np.pi)/hl

print('sque')
hl      = curvedsky.bispec.sque_binnorm(bn,bp,sL)
bl[2,:] = curvedsky.bispec.sque_bin(bn,bp,sL,plm) * np.sqrt(4*np.pi)/hl

print('isos')
hl      = curvedsky.bispec.isos_binnorm(bn,bp)
bl[3,:] = curvedsky.bispec.isos_bin(bn,bp,plm)* np.sqrt(4*np.pi)/hl

# save
np.savetxt('bispec.dat',(np.concatenate((bc[None,:],bl))).T,fmt='%.5e')

