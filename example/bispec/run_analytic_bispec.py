
import basic
import numpy as np
import healpy as hp


import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))


cpmodel = 'modelw'
model = 'RT'
#model = 'GM'
#model = 'SC'

zcmb = 1088.69
zm   = 1.
zs   = [zcmb,zcmb,zcmb]

#zmin, zmax = 0.0001, 1088.69

zmin, zmax = 0.0001, 40.
zn = 50

calc = 'bispec'
#calc = 'bispecsnr'
#calc = 'bispecbin'

btype = 'kkk'
#btype = 'gkk'
#btype = 'ggk'

lmin = 1
lmax = 2048
olmin = lmin
olmax = lmax

bn = 20

D = 'data/'
L = np.linspace(0,olmax,olmax+1)

# Compute matter power spectrum at redshift 0
#k, pk0 = np.loadtxt( D+cpmodel+'/Pk/Pklin.dat', unpack=True )
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.InitPower.set_params(ns=0.965)
pars.set_matter_power(redshifts=[0.], kmax=2.0)

#Linear spectra
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
k, __, pk0 = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
#s8 = np.array(results.get_sigma8())
kn = np.size(k)

z, dz = basic.bispec.zpoints(zmin,zmax,zn)

if btype == 'kkk': 
    dNdz = None
else:
    dNdz = basic.galaxy.dndz_sf(z,2.,1.,zm)

    
# bispectrum
if calc == 'bispec':
    bl0, pb0 = basic.bispec.bispeclens('equi',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    bl1, pb1 = basic.bispec.bispeclens('fold',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    bl2, pb2 = basic.bispec.bispeclens('sque',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    bl3, pb3 = basic.bispec.bispeclens('isos',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    np.savetxt('test_bl_zs'+str(zs[0])+'-'+str(zs[1])+'-'+str(zs[2])+'_'+model+'.dat',np.array((L[1:],bl0,bl1,bl2,bl3,pb0,pb1,pb2,pb3)).T)

    
# binned bispectrum
if calc == 'bispecbin':
    bl, pb = {}, {}
    for shap in ['equi','fold','sque','isos']:
        bc, bl[shap], pb[shap] = basic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zs,lmin,lmax,bn,k,pk0,btype=btype,dNdz=dNdz)
    np.savetxt('test_bl_zs'+str(zs[0])+'-'+str(zs[1])+'-'+str(zs[2])+'_'+model+'_b'+str(bn)+'.dat',np.array((bc,bl['equi'],bl['fold'],bl['sque'],bl['isos'],pb['equi'],pb['fold'],pb['sque'],pb['isos'])).T)

    
# total bispectrum SNR
if calc == 'bispecsnr':

    # load aps
    lmaxs = np.array([100,500,1000,2000,3000])
    #lmaxs = np.array([100,500,1000])
    snr = np.zeros(len(lmaxs))

    for i, lmax in enumerate(lmaxs):
        
        # lensing convergence spectrum
        L = np.linspace(0,lmax,lmax+1)
        clkk = basic.bispec.cl_flat(cpmodel,z,dz,[zmax,zmax],lmax,k,pk0,pktype='T12',cltype='kk')
        
        #clkk = np.zeros(lmax+1)
        #clkk[2:] = np.loadtxt(D+cpmodel+'/cl/fid.dat',unpack=True)[4][:lmax-1]
        
        # noise spectrum if needed
        #nlkk = np.zeros(lmax+1)
        #nlkk[2:] = np.loadtxt(D+'nldd/advact_s6_t1_rlmax4000.dat',unpack=True)[1][:lmax-1]
        #nlkk[2:] = np.loadtxt(D+'nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:lmax-1]
        #clkk[2:] = clkk[2:]/4.*(1.+1./L[2:])**2 + nlkk[2:]#*L[2:]**2*(L[2:]+1.)**2/4.
        
        if btype=='kkk':
            clgg = None
        else:
            clgg = basic.bispec.cl_flat(cpmodel,z,dz,[zmax,zmax],lmax,k,pk0,pktype='T12',cltype='gg',dNdz=dNdz)
            clgg += (np.pi/10800.)**2/30. # shot noise (here assume gal number density per 30 armin squared)

        snr[i] = basic.bispec.bispeclens_snr(cpmodel,model,z,dz,zs,2,lmax,clkk,k,pk0,btype=btype,dNdz=dNdz,cgg=clgg)
        print(lmax,snr[i])

    np.savetxt('snr_'+btype+'_kappa_deproj0_sens2_16000_lT30-3000_lP30-5000_'+str(zn)+'.dat',np.array((lmaxs,snr)).T)


