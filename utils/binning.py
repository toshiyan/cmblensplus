import numpy as np
import sys
#import basic

#////////// Multipole binning //////////

class multipole_binning:

    def __init__(self,n,spc='',lmin=1,lmax=2048,lstart=0):
        '''
        lstart: multipole of cl[:,0], default to 0
        '''
        
        if not isinstance(n,int): 
            sys.exit('n is not integer')

        if not isinstance(spc,str): 
            sys.exit('spc is not str')

        self.n = n
        self.spc = spc
        self.lmin = lmin
        self.lmax = lmax
        self.lstart = lstart
        self.imax = int(self.lmax-self.lstart)  # index of l=lmax
        #self.bp, self.bc = basic.aps.binning(n,[lmin,lmax],spc=spc)
        self.bp, self.bc = binned_ells(n,lmin,lmax,spc)


def binned_ells(n,lmin,lmax,spc=''):
    '''
    return binned multipole edges and centers
    '''

    if lmin<=0 and spc!='':
        sys.exit('ell minimum should be > 0 for spacing')

    if spc == '':
        dl = (lmax-lmin)/n
        bp = np.array( [ lmin+dl*i for i in range(n+1) ] )

    if spc == 'log':
        dl = np.log(lmax/lmin)/n
        bp = np.array( [ lmin*np.exp(dl*i) for i in range(n+1) ] )
    
    if spc == 'log10':
        dl = np.log10(lmax/lmin)/n
        bp = np.array( [ lmin*10.**(dl*i) for i in range(n+1) ] )
    
    if spc == 'p2':
        dl = (np.sqrt(lmax)-np.sqrt(lmin))/n
        bp = np.array( [ (np.sqrt(lmin)+dl*i)**2 for i in range(n+1) ] )
    
    if spc == 'p3':
        dl = (lmax**(1./3.)-lmin**(1./3.))/n
        bp = np.array( [ (lmin**3+dl*i)**3 for i in range(n+1) ] )
    
    bc = np.array( [ (bp[i]+bp[i-1])/2. for i in range(1,n+1)] ) 
    
    return bp, bc
    
        
        
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

    if mb.imax+1 > np.shape(cl)[-1]:
        sys.exit('Multipole number of cl should be equal to or larger than '+str(mb.imax)+', while the input cl has a multipole size of '+str(np.shape(cl)[-1]-1))

    if np.ndim(cl) == 1:
        cb = binning_opt_core(cl[:mb.imax+1],vl[:mb.imax+1],mb)

    if np.ndim(cl) == 2:
        snmax = np.shape(cl)[0]
        cb = np.array([ binning_opt_core(cl[i,:mb.imax+1],vl[:mb.imax+1],mb) for i in range(snmax)])
    
    if np.ndim(cl) == 3:
        snmax = np.shape(cl)[0]
        clnum = np.shape(cl)[1]
        cb = np.array([[ binning_opt_core(cl[i,c,:mb.imax+1],vl[:mb.imax+1],mb) for c in range(clnum)] for i in range(snmax)])

    return cb


def binning_opt_core(cl,vl,mb):

    cb = np.zeros(mb.n)
    
    for i in range(mb.n):
        
        b0, b1 = int(mb.bp[i])-mb.lstart , int(mb.bp[i+1])-mb.lstart
        
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
        
        b0, b1 = int(mb.bp[i])-mb.lstart , int(mb.bp[i+1])-mb.lstart
        
        wl = 1./vl[b0:b1]**2
        
        if np.count_nonzero(wl) > 0:
            wb[i] = 1. / np.sqrt( np.sum(wl[wl!=0]) )
        else:
            wb[i] = 0

    return wb


def binning1(cl,mb):
    """
    dim = 1  ->  cl = [L] and cb = [b]
    dim > 1  ->  cl = [...,L] and cb = [...,b]
    """
    
    if mb.imax+1 > np.shape(cl)[-1]:
        sys.exit('Multipole number of cl should be equal to or larger than '+str(mb.imax)+', while the input cl has a multipole size of '+str(np.shape(cl)[-1]-1))

    vl = np.ones((mb.imax+1))

    if np.ndim(cl) == 1:
        #aps = np.zeros(mb.lmax)
        #aps[mb.lmin:mb.lmax+1] = cl[:mb.imax+1]
        #cb = basic.aps.cl2bcl(mb.n,mb.lmax,aps,lmin=mb.lmin,spc=mb.spc)
        cb = binning_opt_core(cl[:mb.imax+1],vl,mb)

    if np.ndim(cl) == 2:
        snmax = np.shape(cl)[0]
        #aps = np.zeros((snmax,mb.lmax))
        #aps[:,mb.lmin:mb.lmax+1] = cl[:,:mb.imax+1]
        #cb = np.array([basic.aps.cl2bcl(mb.n,mb.lmax,aps[i,:mb.lmax+1],lmin=mb.lmin,spc=mb.spc) for i in range(snmax)])
        cb = np.array([ binning_opt_core(cl[i,:mb.imax+1],vl,mb) for i in range(snmax)])

    if np.ndim(cl) == 3:
        snmax = np.shape(cl)[0]
        clnum = np.shape(cl)[1]
        #aps = np.zeros((snmax,clnum,mb.lmax))
        #aps[:,:,mb.lmin:mb.lmax+1] = cl[:,:,:mb.imax+1]
        #cb = np.array([[basic.aps.cl2bcl(mb.n,mb.lmax,aps[i,c,:mb.lmax+1],lmin=mb.lmin,spc=mb.spc) for c in range(clnum)] for i in range(snmax)])
        cb = np.array([[ binning_opt_core(cl[i,c,:mb.imax+1],vl,mb) for c in range(clnum)] for i in range(snmax)])

    return cb


def binning2(cl,mb0,mb1):

    if mb1.lmin != mb0.lmax+1:
        sys.exit('wrong split')
    if mb1.lmax > np.shape(cl)[-1]-1:
        sys.exit('wrong lmax')

    cb0 = binning1(cl,mb0)
    cb1 = binning1(cl,mb1)
    if np.ndim(cl) == 1:
        return np.concatenate((cb0,cb1))
    if np.ndim(cl) == 2:
        return np.concatenate((cb0,cb1),axis=1)


def binned_spec(mb,fcl,cn=1,doreal=True,opt=False,vl=None,rl=None,lfac=None):
    # for a given array of files, fcl, which containes real (fcl[0]) and sims (fcl[1:]), return realization mean and std of binned spectrum
    # rl is the correction to sim cl
    # lfac is a multipole factor to cl (to be multiplied to both sim and obs cls)
    
    snmax = len(fcl)

    scl = np.array([np.loadtxt(fcl[i],unpack=True)[cn] for i in range(1,snmax)])
    
    if rl is not None:
        scl *= rl[None,:]

    if lfac is not None:
        scl *= lfac[None,:]

    if opt and vl is None:
        Vl = np.std(scl,axis=0)
    else:
        Vl = vl
    
    scb = binning(scl,mb,vl=Vl)

    mcb = np.mean(scb,axis=0)
    vcb = np.std(scb,axis=0)

    if doreal:

        ocl = np.loadtxt(fcl[0],unpack=True)[cn]

        if lfac is not None:
            ocl *= lfac

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






