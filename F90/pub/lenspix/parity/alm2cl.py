
import numpy as np
import healpy as H

simn = 1
lmax = 2000

l  = np.linspace(0,lmax,lmax+1)
cl = np.zeros((3,simn,lmax+1))
for n in range(1,simn+1):
  f = 'alm_r'+str(n)+'.fits'
  tlm = H.fitsfunc.read_alm(f,hdu=1)
  elm = H.fitsfunc.read_alm(f,hdu=2)
  blm = H.fitsfunc.read_alm(f,hdu=3)
  cl[0,n-1,:] = H.sphtfunc.alm2cl(tlm)
  cl[1,n-1,:] = H.sphtfunc.alm2cl(elm)
  cl[2,n-1,:] = H.sphtfunc.alm2cl(blm)

mcl0 = np.mean(cl[0,0:simn],axis=0)
vcl0 = np.std(cl[0,0:simn],axis=0)
mcl1 = np.mean(cl[1,0:simn],axis=0)
vcl1 = np.std(cl[1,0:simn],axis=0)
np.savetxt('mcl.dat',np.transpose((l,mcl0,vcl0,mcl1,vcl1)))

