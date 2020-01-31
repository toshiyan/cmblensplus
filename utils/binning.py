import numpy as np
import sys
import basic

#////////// Multipole binning //////////

class multipole_binning:

    def __init__(self,n,spc='',lmin=1,lmax=2048):
        self.n = n
        self.spc = spc
        self.lmin = lmin
        self.lmax = lmax
        self.bp, self.bc = basic.aps.binning(n,[lmin,lmax],spc=spc)


def binning(cl,mb0,mb1=None):

    if mb1 is None:
        return binning1(cl,mb0)
    else:
        return binning2(cl,mb0,mb1)


def binning1(cl,b):
    """
    dim = 1  ->  cl = (L) and cb = (b)
    dim > 1  ->  cl = (...,L) and cb = (...,b)
    """

    if b.lmax > np.shape(cl)[-1] - 1:
        sys.exit('size of b.lmax is wrong')

    if np.ndim(cl) == 1:
        cb = basic.aps.cl2bcl(b.n,b.lmax,cl[:b.lmax+1],lmin=b.lmin,spc=b.spc)

    if np.ndim(cl) == 2:
        snmax = np.shape(cl)[0]
        cb = np.array([basic.aps.cl2bcl(b.n,b.lmax,cl[i,:b.lmax+1],lmin=b.lmin,spc=b.spc) for i in range(snmax)])

    if np.ndim(cl) == 3:
        snmax = np.shape(cl)[0]
        clnum = np.shape(cl)[1]
        cb = np.array([[basic.aps.cl2bcl(b.n,b.lmax,cl[i,c,:b.lmax+1],lmin=b.lmin,spc=b.spc) for c in range(clnum)] for i in range(snmax)])

    return cb


def binning2(cl,b0,b1):

    if b1.lmin != b0.lmax+1:
        sys.exit('wrong split')
    if b1.lmax > np.shape(cl)[-1]-1:
        sys.exit('wrong lmax')

    cb0 = binning1(cl,b0)
    cb1 = binning1(cl,b1)
    if np.ndim(cl) == 1:
        return np.concatenate((cb0,cb1))
    if np.ndim(cl) == 2:
        return np.concatenate((cb0,cb1),axis=1)



