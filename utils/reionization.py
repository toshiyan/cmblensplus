
import numpy as np, tqdm

# from cmblensplus
import basic

# from cmblensplus/utils
import constants as const
import cosmology

# to avoid scipy constants use
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as spline


def xHlogxH(xe):
    if 1.-xe<1e-6:
        return 0.
    else:
        return (1.-xe)*np.log(1.-xe)


def __Pr__(R,R0,sigma=np.log(2.)):
    return 1./(R*np.sqrt(2*np.pi*sigma**2)) * np.exp(-(np.log(R/R0))**2/(2.*sigma**2))


def __Wth__(x):
    return 3./x**3 * (np.sin(x)-x*np.cos(x))


def __Fk__(k,R0):
    I1 = lambda R: __Pr__(R,R0)*(4*np.pi*R**3/3.)**2*__Wth__(k*R)**2
    I2 = lambda R: __Pr__(R,R0)*4*np.pi*R**3/3.
    neum = quad(I1,0.,2e2)[0]
    deno = quad(I2,0.,2e2)[0]
    return neum/deno


def __Gk_muint__(k,K,Pk):
    I = lambda mu: Pk( np.sqrt(k**2-2*k*K*mu+K**2) )
    return quad(I,-1.,1.)[0]


def __Gk__(k,R0,Pk):
    vFk = np.vectorize(__Fk__)
    vGkint = np.vectorize(__Gk_muint__)
    Ki = np.logspace(-3,1,11)
    dK = Ki[1:]-Ki[:-1]
    I = np.sum(dK*Ki[:-1]**2*vGkint(k,Ki[:-1],Pk)/(2*np.pi)**2*vFk(Ki[:-1],R0))
    return I


def __Ik__(k,R0):
    I1 = lambda R: __Pr__(R,R0) * R**3 * __Wth__(k*R)
    I2 = lambda R: __Pr__(R,R0) * R**3
    neum = quad(I1,0.,2e2)[0]
    deno = quad(I2,0.,2e2)[0]
    return neum/deno


def __cltt__(L,rz,Dz,Hz,Pk,xe,Rz,Kz,bias=6.,zmin=0.,zmax=100.):
    k  = lambda z: L/rz(z)
    Pm = lambda z: Dz(z)**2*Pk(k(z))
    P1 = lambda z: xe(z) * (1-xe(z)) * ( __Fk__(k(z),Rz(z)) + __Gk__(k(z),Rz(z),Pm) )
    P2 = lambda z: ( xHlogxH(xe(z)) * bias * __Ik__(k(z),Rz(z)) - xe(z) )**2 * Pm(z)
    I0 = lambda z: Kz(z)*(P1(z)+P2(z))
    return quad(I0,zmin,zmax)[0]


def xe_sym(z,zre=8.,Dz=4.,f_xe=1.08):
    y = np.power(1.+z,1.5)
    yre = np.power(1.+zre,1.5)
    Dy  = 1.5*np.sqrt(1.+zre)*Dz
    return f_xe*.5*(1.-np.tanh((y-yre)/Dy))


def xe_asym(z,alpha_xe=7.,z_early=20.,zend=6.,f_xe=1.08):
    # Planck 2016 (1605.03507) Eq.(3)
    if z < zend:
        return f_xe
    if z >= zend and z<z_early:
        return f_xe*np.power((z_early-z)/(z_early-zend),alpha_xe)
    if z >= z_early:
        return 0.

    
def optical_depth(xe,H0=70.,Om=.3,Ov=.7,Ob=.0455,w0=-1.,wa=0.,zmin=1e-4,zmax=50,zn=1000):
    
    h0 = H0/100.
    
    # precompute H(z)
    zi  = np.linspace(zmin,zmax,zn)
    Hzi = basic.cosmofuncs.hubble(zi,divc=True,H0=H0,Om=Om,Ov=Ov,w0=w0,wa=wa)
    Hz  = spline( zi, Hzi )

    # define z-integral
    I = lambda z: (1+z)**2/Hz(z) * xe(z)
    
    #Eq.(3.44) of Dodelson's Modern Cosmology and conversion of H(z) unit 1/Mpc -> 1/m
    f = const.sigmaT * (const.rho_c*Ob*h0**2/const.m_p) * const.Mpc2m 
    
    # compute z-integral
    print('optical depth:', quad(I,zmin,zmax)[0] * f) 


def compute_cltt(xe,H0=70.,Om=.3,Ov=.7,Ob=.0455,w0=-1.,wa=0.,ns=.97,As=2e-9,R0=10.,alpha=0.,bias=6.,lmin=1,lmax=3000,ln=100,zmin=1e-4,zmax=50,zn=1000):
    
    cps = {'H0':H0,'Om':Om,'Ov':Ov,'w0':w0,'wa':wa}
    h0  = H0/100.
    
    zi  = np.linspace(zmin,zmax,zn)
    Hzi = basic.cosmofuncs.hubble(zi,divc=True,**cps)
    rzi = basic.cosmofuncs.dist_comoving(zi,**cps)
    Dzi = basic.cosmofuncs.growth_factor(zi,normed=True,**cps)
    Hz  = spline( zi, Hzi )
    rz  = spline( zi, rzi )
    Dz  = spline( zi, Dzi )

    # evolution of bubble size
    if alpha==0.: 
        Rz = lambda z: R0
    else:
        Rz = lambda z: R0*np.power(10.,alpha*(xe(z)-.5))

    # compute linear matter P(k) at z=0
    k, pk0 = cosmology.camb_pk(H0=H0,Om=Om,Ob=Ob,ns=ns,As=As,z=[0.],kmax=20.,minkh=1e-4,maxkh=10,npoints=500)
    Pk = spline(k,pk0)
    
    # Kz
    Kz = lambda z: (const.sigmaT*(const.rho_c*Ob*h0**2/const.m_p)*const.Mpc2m)**2 * (1+z)**4/rz(z)**2/Hz(z)

    # compute cltt
    l  = np.linspace(lmin,lmax,ln)
    cl = np.zeros(ln)
    for i, L in enumerate(tqdm.tqdm(l)):
        cl[i] = __cltt__(L,rz,Dz,Hz,Pk,xe,Rz,Kz,bias=bias)

    return l, cl

