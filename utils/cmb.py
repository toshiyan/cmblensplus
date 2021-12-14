# external

import numpy as np
import healpy as hp
import sys
import tqdm

# from cmblensplus/wrap
import basic
import curvedsky

# from cmblensplus/utils
import constant as c
import misctools

# for pickle (to be removed)
if sys.version_info[:3] > (3,0):
    import pickle
elif sys.version_info[:3] > (2,5,2):
    import cPickle as pickle


#////////// Constants //////////#

# CMB temperature in K
Tcmb = c.Tcmb

# CMB types
mtype = ['T','E','B']


#////////// Frequency Spectrum //////////#

def Int_Planck(nu,T): 
    '''
    Black body I(nu)
    nu [GHz]
    '''
    return 2*c.h*(nu*1e9)**3/c.c**2 * (np.exp(c.h*(nu*1e9)/c.k/T)-1.)**(-1)


def Int_dust(nu,Td=19.6,beta=1.53):
    '''
    Dust I(nu) as a modified Black body
    Td, beta: values from arXiv 1801.04945
    nu [GHz]
    '''
    return (nu*1e9)**beta*Int_Planck(nu,Td)


def fnu_dust(nu,Td=34.,beta=2.,nu_c=4955.,alpha=2.,bc=1e-64):
    '''
    f(nu) defined in Eq.(D.4) of arXiv:1502.01591
    nu [GHz]
    '''
    fnu = np.zeros(len(nu))
    for i, nu_i in enumerate(nu):
        if nu_i <= nu_c:
            fnu[i] = (nu_i*1e9)**beta*Int_Planck(nu_i,Td)*c.c**2/(2*c.h)
        else:
            fnu[i] = (nu_c*1e9)**beta*Int_Planck(nu_c,Td)*c.c**2/(2*c.h)*(nu_i/nu_c)**(-alpha)
    return fnu*bc
    

def RJfunc(nu,T):
    """
    nu [GHz]
    """
    return 2*(nu*1e9)**2/c.c**2 * c.k*T


def Int_Planck_deriv(nu,T):
    """
    nu [GHz]
    """
    nu0 = nu*1e9
    return Int_Planck(nu,T) * (np.exp(c.h*nu0/c.k/T)*c.h*nu0/c.k/T**2)/(np.exp(c.h*nu0/c.k/T)-1.)


def T2y(nu):
    """
    g(nu) factor to convert from T to y
    nu [GHz]
    """
    x = c.h*(nu*1e9)/(c.k*c.Tcmb*1e-6)
    return x*np.cosh(x/2.)/np.sinh(x/2.) - 4.


def unit_conversion(freq,conv): #https://lambda.gsfc.nasa.gov/data/planck/early_sources/explanatory_supplement.pdf Table 4

    unit = {}
    unit['30']  = {'center':28.5, 'Kcmb->KRJ':0.979328,     'KRJ->MJy/sr':24.845597}
    unit['44']  = {'center':44.1, 'Kcmb->KRJ':0.95121302,   'KRJ->MJy/sr':59.666236}
    unit['70']  = {'center':70.3, 'Kcmb->KRJ':0.88140690,   'KRJ->MJy/sr':151.73238}
    unit['100'] = {'center':100,  'Kcmb->KRJ':0.76581996,   'KRJ->MJy/sr':306.81118}
    unit['143'] = {'center':143,  'Kcmb->KRJ':0.59714682,   'KRJ->MJy/sr':627.39818}
    unit['217'] = {'center':217,  'Kcmb->KRJ':0.31573332,   'KRJ->MJy/sr':1444.7432}
    unit['353'] = {'center':353,  'Kcmb->KRJ':0.071041398,  'KRJ->MJy/sr':3823.1434}
    unit['545'] = {'center':545,  'Kcmb->KRJ':0.0059757149, 'KRJ->MJy/sr':9113.0590}
    unit['857'] = {'cenrer':857,  'Kcmb->KRJ':9.6589431e-05,'KRJ->MJy/sr':22533.716}

    if conv in ['Kcmb->KRJ','KRJ->MJy/sr']:
        return unit[freq][conv]
    elif conv == 'Kcmb->MJy/sr':
        return unit[freq]['Kcmb->KRJ']*unit[freq]['KRJ->MJy/sr']


#////////// Gaussian beam //////////#

def beam(theta,lmax,inv=True):
    # inverse of beam function
    L = np.linspace(0,lmax,lmax+1)
    if inv:
        return np.exp(L*(L+1)*(theta*c.ac2rad)**2/(16.*np.log(2.)))
    else:
        return np.exp(-L*(L+1)*(theta*c.ac2rad)**2/(16.*np.log(2.)))


#////////// Noise spectrum //////////#

def nl_white(sigma,theta,lmax):

    return (sigma*c.ac2rad/c.Tcmb)**2*beam(theta,lmax)**2


#////////// Angular power spectrum /////////#

def aps(rlz,lmax,falm,odd=True,w2=1.,mtype=['T','E','B'],fname=None,skip_rlz=[],loadcls=True,verbose=True,overwrite=False):
    '''
    Compute CMB aps (TT,EE,BB,TE,TB,EB)
    '''

    if odd:
        cn = 6
    else:
        cn = 4

    cls = np.zeros((len(rlz),cn,lmax+1))

    for ii, i in enumerate(tqdm.tqdm(rlz,ncols=100,desc='aps:')):

        if i in skip_rlz: continue

        if fname is not None and misctools.check_path(fname[i],verbose=verbose,overwrite=overwrite):

            if loadcls:
                cls[ii,:,:] = np.loadtxt(fname[i],unpack=True)[1:,:]

        else:
        
            # compute cls
            if 'T' in mtype:  
                Talm = pickle.load(open(falm['T'][i],"rb"))
                cls[ii,0,:] = curvedsky.utils.alm2cl(lmax,Talm)

            if 'E' in mtype:  
                Ealm = pickle.load(open(falm['E'][i],"rb"))
                cls[ii,1,:] = curvedsky.utils.alm2cl(lmax,Ealm)

            if 'B' in mtype:  
                Balm = pickle.load(open(falm['B'][i],"rb"))
                cls[ii,2,:] = curvedsky.utils.alm2cl(lmax,Balm)

            if 'T' in mtype and 'E' in mtype:  
                cls[ii,3,:] = curvedsky.utils.alm2cl(lmax,Talm,Ealm)

            if odd and 'T' in mtype and 'B' in mtype:
                cls[ii,4,:] = curvedsky.utils.alm2cl(lmax,Talm,Balm)

            if odd and 'E' in mtype and 'B' in mtype:
                cls[ii,5,:] = curvedsky.utils.alm2cl(lmax,Ealm,Balm)

            # correct normalization
            if w2 != 1.:
                cls[ii,:,:] /= w2

        # save to file
        if fname is not None and overwrite:
            L = np.linspace(0,lmax,lmax+1)
            np.savetxt(fname[i],np.concatenate((L[None,:],cls[ii,:,:])).T)

    return cls


def apsx(rlz,lmax,falm,galm,verbose=True,overwrite=False,mtype=['T','E','B']):
    '''
    Compute CMB aps (T0T1,E0E1,B0B1)
    '''

    cls = np.zeros((len(rlz),3,lmax+1))

    for ii, i in enumerate(tqdm.tqdm(rlz,ncols=100,desc='apsx:')):

        if 'T' in mtype:  
            Talm = pickle.load(open(falm['T'][i],"rb"))[:lmax+1,:lmax+1]
            talm = pickle.load(open(galm['T'][i],"rb"))[:lmax+1,:lmax+1]
            cls[ii,0,:] = curvedsky.utils.alm2cl(lmax,Talm,talm)

        if 'E' in mtype:  
            Ealm = pickle.load(open(falm['E'][i],"rb"))[:lmax+1,:lmax+1]
            ealm = pickle.load(open(galm['E'][i],"rb"))[:lmax+1,:lmax+1]
            cls[ii,1,:] = curvedsky.utils.alm2cl(lmax,Ealm,ealm)

        if 'B' in mtype:  
            Balm = pickle.load(open(falm['B'][i],"rb"))[:lmax+1,:lmax+1]
            balm = pickle.load(open(galm['B'][i],"rb"))[:lmax+1,:lmax+1]
            cls[ii,2,:] = curvedsky.utils.alm2cl(lmax,Balm,balm)

    return cls



def read_camb_cls(fname,ftype='scal',output='',skiprows=1,unpack=True):
    
    if ftype == 'scal': # unlensed scal Cls
        ll, TT, EE, TE, PP, TP = np.loadtxt(fname,skiprows=skiprows,unpack=unpack)
    
    elif ftype == 'lens': # lensed cls
        ll, TT, EE, BB, TE = np.loadtxt(fname,skiprows=skiprows,unpack=unpack)

    s = ll*(ll+1.)*c.Tcmb**2/(2*np.pi)
    TT /= s
    EE /= s
    TE /= s
    if ftype == 'scal':
        PP /= ll**4*c.Tcmb**2
        TP /= ll**3*c.Tcmb**2
    elif ftype == 'lens':
        BB /= s

    TT = np.insert(TT,0,np.array([0.,0.]))
    EE = np.insert(EE,0,np.array([0.,0.]))
    TE = np.insert(TE,0,np.array([0.,0.]))
    if ftype == 'scal':
        PP = np.insert(PP,0,np.array([0.,0.]))
        TP = np.insert(TP,0,np.array([0.,0.]))
    elif ftype == 'lens':
        BB = np.insert(BB,0,np.array([0.,0.]))

    if ftype == 'scal':
        if output=='array':
            return np.array((TT,EE,TE,PP,TP))
        else:
            return TT, EE, TE, PP, TP
    elif ftype == 'lens':
        if output=='array':
            return np.array((TT,EE,BB,TE))
        else:
            return TT, EE, BB, TE


#////////// Window function /////////#

def wfactor(w,nmax=5):
    
    # normalization correction due to window (n=0 will be replaced below)
    wn = np.array( [ np.average(w**n) for n in range(nmax) ] )
    
    # binary mask for n=0
    if not np.isscalar(w):
        m = w.copy()
        m[m!=0.] = 1.
        wn[0] = np.average(m)
    
    return wn

#////////// Map -> Alm //////////#
def map2alm_curvedsky(lmax,fmap,mtype=['T','E','B'],w=1.,bl=None,uKcmb=False):

    if uKcmb:
        norm = c.Tcmb
    else:
        norm = 1. 
        
    # load map (in unit of Tcmb for default)
    if 'T' in mtype:
        Tmap = w * hp.fitsfunc.read_map(fmap,field=0,verbose=False)/norm
    if 'E' in mtype or 'B' in mtype:
        Qmap = w * hp.fitsfunc.read_map(fmap,field=1,verbose=False)/norm
        Umap = w * hp.fitsfunc.read_map(fmap,field=2,verbose=False)/norm
    
    # get nside
    nside = hp.pixelfunc.get_nside(Tmap)

    # map to alm
    alm = {}
    if 'T' in mtype:
        alm['T'] = curvedsky.utils.hp_map2alm(nside,lmax,lmax,Tmap)
    if 'E' in mtype or 'B' in mtype:
        alm['E'], alm['B'] = curvedsky.utils.hp_map2alm_spin(nside,lmax,lmax,2,Qmap,Umap)

    # beam deconvolution
    if bl is not None:
        for m in mtype:
            alm[m] /= bl[:,None]

    return alm

