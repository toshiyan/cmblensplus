import numpy as np
import sys
import basic

#////////// Multipole binning //////////

class multipole_binning:

    def __init__(self,n,spc='',lmin=1,lmax=2048):
        
        if not isinstance(n,int): 
            sys.exit('n is not integer')

        if not isinstance(spc,str): 
            sys.exit('spc is not str')

        self.n = n
        self.spc = spc
        self.lmin = lmin
        self.lmax = lmax
        self.bp, self.bc = basic.aps.binning(n,[lmin,lmax],spc=spc)


def binning(cl,mb0,mb1=None,vl=None):

    if vl is None:
        if mb1 is None:
            return binning1(cl,mb0)
        else:
            return binning2(cl,mb0,mb1)
    else:
        # binning of power spectrum with un-equal weights from variance
        return binning_opt(cl,vl,mb0)


def binning_opt(cl,vl,mb):

    if mb.lmax > np.shape(cl)[-1] - 1:
        sys.exit('size of mb.lmax is wrong: '+str(mb.lmax)+', '+str(np.shape(cl)[-1]-1))

    if np.ndim(cl) == 1:
        cb = binning_opt_core(cl[:mb.lmax+1],vl[:mb.lmax+1],mb)

    if np.ndim(cl) == 2:
        snmax = np.shape(cl)[0]
        cb = np.array([ binning_opt_core(cl[i,:mb.lmax+1],vl[:mb.lmax+1],mb) for i in range(snmax)])
    
    if np.ndim(cl) == 3:
        snmax = np.shape(cl)[0]
        clnum = np.shape(cl)[1]
        cb = np.array([[ binning_opt_core(cl[i,c,:mb.lmax+1],vl[:mb.lmax+1],mb) for c in range(clnum)] for i in range(snmax)])

    return cb


def binning_opt_core(cl,vl,mb):

    cb = np.zeros(mb.n)
    
    for i in range(mb.n):
        
        b0, b1 = int(mb.bp[i]) , int(mb.bp[i+1])
        
        cs = cl[b0:b1]
        wl = 1./vl[b0:b1]**2
        
        if np.count_nonzero(wl) > 0:
            cb[i] = np.sum(wl[wl!=0]*cs[wl!=0]) / np.sum(wl[wl!=0])
        else:
            cb[i] = 0

    return cb


def binning_opt_weight(vl,mb):

    wb = np.zeros(mb.n)
    
    for i in range(mb.n):
        
        b0, b1 = int(mb.bp[i]) , int(mb.bp[i+1])
        
        wl = 1./vl[b0:b1]**2
        
        if np.count_nonzero(wl) > 0:
            wb[i] = 1. / np.sqrt( np.sum(wl[wl!=0]) )
        else:
            wb[i] = 0

    return wb


def binning1(cl,b):
    """
    dim = 1  ->  cl = [L] and cb = [b]
    dim > 1  ->  cl = [...,L] and cb = [...,b]
    """

    if b.lmax > np.shape(cl)[-1] - 1:
        sys.exit('size of b.lmax is wrong: '+str(b.lmax)+', '+str(np.shape(cl)[-1]-1))

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


def binned_spec(mb,fcl,cn=1,doreal=True,opt=False,vl=None,rl=None):
    # for a given array of files, fcl, which containes real (fcl[0]) and sims (fcl[1:]), return realization mean and std of binned spectrum
    # rl is the correction to cl
    
    snmax = len(fcl)

    scl = np.array([np.loadtxt(fcl[i],unpack=True)[cn] for i in range(1,snmax)])
    
    if rl is not None:
        scl *= rl[None,:]

    if opt and vl is None:
        Vl = np.std(scl,axis=0)
    else:
        Vl = vl
    
    scb = binning(scl,mb,vl=Vl)

    mcb = np.mean(scb,axis=0)
    vcb = np.std(scb,axis=0)

    if doreal:

        ocl = np.loadtxt(fcl[0],unpack=True)[cn]
        ocb = binning(ocl,mb,vl=Vl)

        if opt and vl is None:
            return mcb, vcb, scb, ocb, Vl
        else:
            return mcb, vcb, scb, ocb

    else:
   
        if opt and vl is None:
            return mcb, vcb, scb, Vl
        else:
            return mcb, vcb, scb






