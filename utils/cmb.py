import numpy as np
import sys
import constants as c
import basic
import curvedsky

if sys.version_info[:3] > (3,0):
  import pickle
elif sys.version_info[:3] > (2,5,2):
  import cPickle as pickle


#////////// Frequency Spectrum //////////#

def Int_Planck(nu,T): 
    '''
    Black body I(nu)
    nu [GHz]
    '''
    nu0 = nu*1e9
    return 2*c.h*nu0**3/c.c**2 * (np.exp(c.h*nu0/c.kB/T)-1.)**(-1)


def Int_dust(nu,Td=19.6,beta=1.53): 
    '''
    Dust I(nu) as a modified Black body
    Td, beta: values from arXiv 1801.04945
    nu [GHz]
    '''
    return (nu*1e9)**beta*Int_Planck(nu,Td)


def RJfunc(nu,T):
    """
    nu [GHz]
    """
    return 2*(nu*1e9)**2/c.c**2 * c.kB*T


def Int_Planck_deriv(nu,T):
    """
    nu [GHz]
    """
    nu0 = nu*1e9
    return Int_Planck(nu,T) * (np.exp(c.h*nu0/c.kB/T)*c.h*nu0/c.kB/T**2)/(np.exp(c.h*nu0/c.kB/T)-1.)


#////////// Gaussian beam //////////#

def beam(theta,lmax):

    ac2rad = np.pi/10800.
    L      = np.linspace(0,lmax,lmax+1)
    beam   = np.exp(L*(L+1)*(theta*ac2rad)**2/(16.*np.log(2.)))
    return beam


#////////// Noise spectrum //////////#

def nl_white(sigma,theta,lmax,Tcmb=2.72e6):

    return (sigma*ac2rad/Tcmb)**2*beam(theta,lmax)**2


#////////// Angular power spectrum /////////#

def aps(snmin,snmax,lmax,falm,verbose=True):
    '''
    Compute CMB aps (TT,EE,BB,TE)
    '''

    if verbose: print('aps')

    sn  = snmax - snmin + 1
    #cbs = np.zeros((sn,4,bn))
    cls = np.zeros((sn,4,lmax+1))

    for i in range(snmin,snmax+1):

        if verbose: print(i, end=' ')
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

        #for j in range(4):
        #    cbs[ii,j,:] = basic.aps.cl2bcl(bn,lmax,cls[ii,j,:],spc=binspc)

    return cls#, cbs



