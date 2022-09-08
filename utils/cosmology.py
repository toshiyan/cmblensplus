
# general
import numpy as np
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.interpolate import interp2d
from astropy import units as u
import hankel
import camb

# cmblensplus/wrap/
import basic

# cmblensplus/utils/
import constant
import cmb
import analysis as ana


#////////// Constants //////////#

# Comoving weight function for computing angular power spectrum of large-scale structure mass distribution

def window_cib(chi,z,nu,z_c=2.,sigma_z=2.,bc=.528654,dim='',weight_only=False):
    '''
    The window has a unit of 1/Mpc (see Eq.(D.3) of arXiv:1502.01591). To pass into camb python as a source count, 1/H(z) should be further multiplied. 
    rz: comoving distance [Mpc]
    z : redshift
    nu: frequency [GHz]
    '''
    if weight_only:
        w_cib = bc * (chi**2/(1.+z)**2) * cmb.fnu_dust(nu*(1+z),Td=34.,beta=2.) 
    else:
        w_cib = bc * (chi**2/(1.+z)**2) * np.exp(-(z-z_c)**2/(2*sigma_z**2)) * cmb.fnu_dust(nu*(1+z),Td=34.,beta=2.) 

    if dim == '':
        w_cib *= constant.Jysr2uK(nu)/constant.Tcmb
    elif dim == 'uK':
        w_cib *= constant.Jysr2uK(nu)
    elif dim == 'Jysr': 
        w_cib *= 1.
        

    return w_cib


#def window_gal(Hz,chi,fNz):
#    '''
#    # fNz is a normalized z-distribution so that its integral over z becomes a fractional number of galaxies at the bin
#    # For a single bin, the integral of fNz over z becomes unity. 
#    The window has a unit of 1/Mpc
#    '''
#    return (Hz/chi) * fNz


def window_kap(z,zcmb=1088.69,**cps):
    '''
    The window has a unit of 1/Mpc (see Eq.(D.3) of arXiv:1502.01591). To pass into camb python as a source count, 1/H(z) should be further multiplied.
    '''
    chi  = basic.cosmofuncs.dist_comoving(z,**cps)
    chis = basic.cosmofuncs.dist_comoving([zcmb],**cps)[0]
    return 1.5 * cps['Om']*(cps['H0']/constant.C)**2 * (1.+z) * chi*(chis-chi)/chis


# //// CAMB Cl //// #

def calc_cmb_aps(lmax,cltype=[],H0=67.5,ombh2=0.022,omch2=0.122,mnu=0.06,omk=0,tau=0.06,As=2e-9,ns=0.965,r=0,lens_ac=2):
    
    #Set up a new set of parameters for CAMB
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk, tau=tau)
    pars.InitPower.set_params(As=As, ns=ns, r=r)
    pars.set_for_lmax(lmax, lens_potential_accuracy=lens_ac);
    
    #calculate results for these parameters
    results = camb.get_results(pars)

    #read lensed CMB aps (TT,EE,BB,TE)
    ucl, lcl = None, None
    l = np.linspace(0,lmax,lmax+1)
    if 'unlens' in cltype:
        ucl = results.get_unlensed_scalar_cls().T[:,:lmax+1] * 2*np.pi/(l[None,:]**2+l[None,:]+1e-30)
    if 'lens' in cltype:
        lcl = results.get_lensed_scalar_cls().T[:,:lmax+1] * 2*np.pi/(l[None,:]**2+l[None,:]+1e-30)
    return ucl, lcl
        

def calc_phi_aps(Lmax,H0=67.5,ombh2=0.022,omch2=0.122,mnu=0.06,omk=0,As=2e-9,ns=0.965,ac=1,lsamp=1,lens_ac=5):

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_for_lmax(Lmax, lens_potential_accuracy=lens_ac)
    pars.NonLinear = camb.model.NonLinear_both
    pars.Accuracy.AccuracyBoost = ac
    pars.Accuracy.lSampleBoost = lsamp
    
    results  = camb.get_results(pars)
    
    return results.get_lens_potential_cls(Lmax)


def calc_lss_aps(Lmax,zi,w,bz=1.,H0=67.5,ombh2=0.022,omch2=0.122,mnu=0.06,omk=0,As=2e-9,ns=0.965,ac=1,lsamp=1,lens_ac=5):

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=omk)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_for_lmax(Lmax, lens_potential_accuracy=lens_ac)

    pars.Want_CMB = False
    pars.NonLinear = camb.model.NonLinear_both
    pars.Accuracy.AccuracyBoost = ac
    pars.Accuracy.lSampleBoost = lsamp
    
    #//// Set up W(z) window functions ////#
    klist = list(w)
    
    if 'cib' in klist:
        tracers = [ camb.sources.SplinedSourceWindow( z=zi, W=w['cib'], dlog10Ndm=.4, bias=np.sum(w['cib']*(zi[1]-zi[0])) ) ]
    
    for m in klist: # add galaxies
        if m == 'cib':  continue
        tracers += [ camb.sources.SplinedSourceWindow( z=zi, W=w[m], dlog10Ndm=0, bias_z=bz ) ]
    pars.SourceWindows = tracers

    # turning off GR corrections as they impact on large-scale (ell<20) for galaxy and should not present in CIB
    pars.SourceTerms.counts_redshift = False 
    pars.SourceTerms.counts_velocity = False
    pars.SourceTerms.counts_timedelay = False
    pars.SourceTerms.counts_ISW = False
    pars.SourceTerms.counts_potential = False
    camb_list = np.concatenate((np.array(['P']),np.array(list([ 'W'+str(i) for i in range(1,len(w)+1) ]))))
        
    results   = camb.get_results(pars)
    
    return camb_list, results.get_source_cls_dict()

    

# //// Limber Cl //// #

def cl_limber(L,z,dz,W0,W1,k,pk0,ln=200,**cps):
    '''
    The weights, W0 and W1, have a unit of 1/Mpc. 
    k:   Wavenumber in unit of h/Mpc
    pk0: P(k,z=0) in unit of (Mpc/h)^3
    '''
    
    h0  = cps['H0']*1e-2
    Pk  = spline(k*h0,pk0/h0**3)

    Hzi = basic.cosmofuncs.hubble(z,divc=True,**cps)
    chi = basic.cosmofuncs.dist_comoving(z,**cps)
    Dzi = basic.cosmofuncs.growth_factor(z,normed=True,**cps)

    ll  = np.linspace(L[0],L[-1],ln)
    cli = np.zeros(ln)

    for i, li in enumerate(ll):
    
        integ  = (dz/Hzi) * (Hzi/chi)**2 * (W0/Hzi)*(W1/Hzi) * Dzi**2*Pk((li+.5)/chi)  
        cli[i] = np.sum(integ[:-1])
    
    return  spline ( ll, cli )(L)


# Matter P(k)

def camb_pk(H0=70.,Om=.3,Ob=0.0455,ns=.96,As=2e-9,z=[0.],kmax=20.,minkh=1e-4,maxkh=10,npoints=500):
    # compute linear matter P(k) with CAMB
    h0 = H0*.01
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=Ob*h0**2, omch2=(Om-Ob)*h0**2)
    pars.InitPower.set_params(ns=ns,As=As)
    pars.set_matter_power(redshifts=z, kmax=kmax)
    results = camb.get_results(pars)
    k, __, pk0 = results.get_matter_power_spectrum(minkh=minkh, maxkh=maxkh, npoints=npoints)
    return  k*h0, pk0/h0**3


def Pmk(k,z,k0=0.05,H0=70.,ns=.96,Om=.27,As=2e-9):
    fc = 2*np.pi**2/np.power(H0,3+ns)
    return (8.*np.pi**2/(25.*Om*(H0*.01)**2))*As*k*np.power(k/k0,ns-1.)*Tk_BBKS(k,Om=Om,H0=H0)**2*basic.cosmofuncs.growth_factor(z,H0=H0,Om=Om)**2


def Tk_BBKS(k,Om=.27,H0=70.):
    # Transfer function (BBKS)
    q  = k/(Om*(H0*.01)**2)
    return np.log(1.+2.34*q)/(2.34*q) / np.power(1.+3.89*q+(16.2*q)**2+(5.47*q)**3+(6.71*q)**4,.25)


# Correlation function

def pk2corrfunc(k,pk,r,h=0.001):
    # compute correlation function for a given matter power spectrum
    #r = np.arange(rmin,rmax,dr)
    # power spectrum
    pkfunc = spline(k, k**0.5/(2**1.5*np.pi**1.5)*pk )
    # setup hankel transform
    ht = hankel.HankelTransform(nu=0.5,h=h)
    # compute correlation function
    XI = ht.transform(pkfunc, r, False, inverse=True) / np.sqrt(r)
    return spline(r,XI)


def IntegPJ0(k,pk0,h=.001,verbose=False):
    # construct a function: F(r) = int[kP(k)*J_0(kr)/2pi]
    r = basic.cosmofuncs.dist_comoving( np.linspace(1e-4,5.,2000) ) * .05
    if verbose: print('compute using r between',min(r),max(r),r[1]-r[0])
    # setup hankel transform
    ht = hankel.HankelTransform(nu=0.,h=h)
    # not use Hankel transform for r = 0 due to numerical accuracy
    PJ0 = integrate.quad( spline(k, k*pk0) , 1e-4, 4. )[0]/(2*np.pi)
    # for r /= 0
    PJR = ht.transform( spline(k, pk0) , r, False, inverse=True)/(2*np.pi)
    # replace r=0 with the exact result
    PJR[0] = PJ0
    return spline(r,PJR)


def IntegPJ(k,zs,pk,h=.001):
    # construct a function: F(r) = int[kP(k)*J_0(kr)/2pi]
    r = basic.cosmofuncs.dist_comoving( np.linspace(1e-4,5.,2000) ) * .05
    # setup hankel transform
    ht = hankel.HankelTransform(nu=0.,h=h)
    PJR = np.zeros((len(zs),len(r)))
    for zi, z in enumerate(zs):
        # for r /= 0
        PJR[zi,:] = ht.transform( spline(k, pk[zi]) , r, False, inverse=True)/(2*np.pi)
        # replace r=0 with the exact result
        PJR[zi,0] = integrate.quad( spline(k, k*pk[zi]) , 1e-4, 4. )[0]/(2*np.pi)
    return interp2d(r,zs,PJR,kind='cubic')


def phi2deltam_coeff(z,**cp):
    Dz = basic.cosmofuncs.growth_factor([z],normed=True,**cp)[0]
    return  1.5*cp['Om']*(cp['H0']/3e5)**2*(1.+z)*Dz(z)


#//// Forecast ////#

def Fisher_Matrix(L,dCdp=None,iC=None,dlnCdp=None,fsky=1):
    return ana.Fisher_Matrix(L,dCdp=dCdp,iC=iC,dlnCdp=dlnCdp,fsky=fsky)


def Fisher_2Dcontour(F,i=0,j=1,display=False):
    return ana.Fisher_2Dcontour(F,i=i,j=j,display=display)
    
