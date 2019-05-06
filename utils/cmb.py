import numpy as np
import sys
import basic
import curvedsky
from memory_profiler import profile

if sys.version_info[:3] > (3,0):
  import pickle
elif sys.version_info[:3] > (2,5,2):
  import cPickle as pickle


@profile
def aps(snmax,bn,binspc,lmax,falm):
  '''
  Compute CMB aps (TT,EE,BB,TE)
  '''

  cbs  = np.zeros((snmax,4,bn))
  cls  = np.zeros((snmax,4,lmax+1))

  for i in range(snmax):

    print('load alm', i)

    # load cmb alms
    Talm = pickle.load(open(falm['T'][i],"rb"))
    Ealm = pickle.load(open(falm['E'][i],"rb"))
    Balm = pickle.load(open(falm['B'][i],"rb"))

    # compute cls
    cls[i,0,:] = curvedsky.utils.alm2cl(lmax,Talm)
    cls[i,1,:] = curvedsky.utils.alm2cl(lmax,Ealm)
    cls[i,2,:] = curvedsky.utils.alm2cl(lmax,Balm)
    cls[i,3,:] = curvedsky.utils.alm2cl(lmax,Talm,Ealm)

    for j in range(4):
      cbs[i,j,:] = basic.aps.cl2bcl(bn,lmax,cls[i,j,:],spc=binspc)

  return cls, cbs

