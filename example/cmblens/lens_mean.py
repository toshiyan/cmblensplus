# example of mean-fiels bias calculation
import numpy as np
import curvedsky
import pickle

simn = 100
lmax = 3000

# need reconstructed lensing estimator
fname = ['test.dat']
fmean = ['mean.dat']

for i in range(1):

  cl = np.zeros((2,lmax+1))

  mfg = 0
  mfc = 0

  for I in range(100):

    if I == i:  continue
    print(I)

    glm, clm = pickle.load(open(fname[I],"rb"))
    mfg += glm/np.float(simn-1)
    mfc += clm/np.float(simn-1)

  # compute mf cls
  cl[0,:] = curvedsky.utils.alm2cl(lmax,mfg)
  cl[1,:] = curvedsky.utils.alm2cl(lmax,mfc)

  np.savetxt(fmean[i],np.concatenate((eL[None,:],cl)).T)

