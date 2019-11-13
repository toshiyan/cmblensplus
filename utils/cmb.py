import numpy as np
import sys
import basic
import curvedsky

if sys.version_info[:3] > (3,0):
  import pickle
elif sys.version_info[:3] > (2,5,2):
  import cPickle as pickle


def aps(snmin,snmax,bn,binspc,lmax,falm):
    '''
    Compute CMB aps (TT,EE,BB,TE)
    '''

    sn  = snmax - snmin
    cbs = np.zeros((sn,4,bn))
    cls = np.zeros((sn,4,lmax+1))

    for i in range(snmin,snmax):

        print('load alm', i)
        ii = i - snmin

        # load cmb alms
        Talm = pickle.load(open(falm['T'][i],"rb"))
        Ealm = pickle.load(open(falm['E'][i],"rb"))
        Balm = pickle.load(open(falm['B'][i],"rb"))

        # compute cls
        cls[ii,0,:] = curvedsky.utils.alm2cl(lmax,Talm)
        cls[ii,1,:] = curvedsky.utils.alm2cl(lmax,Ealm)
        cls[ii,2,:] = curvedsky.utils.alm2cl(lmax,Balm)
        cls[ii,3,:] = curvedsky.utils.alm2cl(lmax,Talm,Ealm)

        for j in range(4):
            cbs[ii,j,:] = basic.aps.cl2bcl(bn,lmax,cls[ii,j,:],spc=binspc)

    return cls, cbs


# Noise
def nl_white(sigma,theta,lmax,Tcmb=2.72e6):
    ac2rad = np.pi/10800.
    noise  = (sigma*ac2rad/Tcmb)**2
    L      = np.linspace(0,lmax,lmax+1)
    beam   = np.exp(L*(L+1)*(theta*ac2rad)**2/(8.*np.log(2.)))
    return noise*beam


# beam
def beam(theta,lmax):
    ac2rad = np.pi/10800.
    L      = np.linspace(0,lmax,lmax+1)
    beam   = np.exp(-L*(L+1)*(theta*ac2rad)**2/(16.*np.log(2.)))
    #Lmax   = int(10./(theta*ac2rad))
    #beam[Lmax:] = beam[min(lmax,Lmax-1)]
    return beam


