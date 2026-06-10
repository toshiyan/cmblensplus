import numpy as np
import ducc0
from .. import libcurvedsky


def alm2map(nside: int, alm: Array, nthreads=0) -> Array:

    # extract maximum multipole
    lmax = len(alm[:,0]) - 1
    mmax = len(alm[0,:]) - 1

    hb = ducc0.healpix.Healpix_Base(nside, "RING")

    # convert alm 2D array to alm 1D array
    lmpy  = int((lmax+1)*(lmax+2)/2.)
    almpy = libcurvedsky.utils.lm_healpix2healpy(lmax,alm,lmpy)

    # preparation for SHT
    alm_ducc = np.asarray(almpy, dtype=np.complex128).reshape(1, -1)
    mp = np.empty((1, hb.npix()), dtype=np.float64)

    sht_info = hb.sht_info()

    # SHT
    ducc0.sht.synthesis(alm=alm_ducc,map=mp,spin=0,lmax=lmax,mmax=mmax,nthreads=nthreads,**sht_info)
    
    return np.asarray(mp[0], dtype=np.float64)


def alm2map_spin(nside: int, spin: int, alm: Array, nthreads=0) -> Array:
    # expected input:
    # alm[0,l,m] = first spin component, e.g. E
    # alm[1,l,m] = second spin component, e.g. B

    # extract maximum multipole
    lmax = len(alm[0,:,0]) - 1
    mmax = len(alm[0,0,:]) - 1

    hb  = ducc0.healpix.Healpix_Base(nside, "RING")

    # convert each 2D alm[l,m] to healpy packed 1D alm
    lmpy   = int((lmax+1)*(lmax+2)/2.)
    almpy0 = libcurvedsky.utils.lm_healpix2healpy(lmax, alm[0], lmpy)
    almpy1 = libcurvedsky.utils.lm_healpix2healpy(lmax, alm[1], lmpy)

    # preparation for ducc0 SHT
    alm_ducc = np.asarray([almpy0, almpy1], dtype=np.complex128)
    mp = np.empty((2, hb.npix()), dtype=np.float64)
    
    sht_info = hb.sht_info()

    # SHT
    ducc0.sht.synthesis(alm=alm_ducc,map=mp,spin=spin,lmax=lmax,mmax=mmax,nthreads=nthreads,**sht_info)

    return np.asarray(mp, dtype=np.float64)

    
def map2alm(lmax: int, mmax: int, map_in: Array, nthreads=0, maxiter=0) -> Array:

    nside = int(np.sqrt(len(map_in) / 12))

    hb = ducc0.healpix.Healpix_Base(nside, "RING")
    sht_info = hb.sht_info()
    
    ducc_map = np.asarray(map_in, dtype=np.float64).reshape(1, -1)
    
    almpy = ducc0.sht.adjoint_synthesis(map=ducc_map,spin=0,lmax=lmax,mmax=mmax,nthreads=nthreads,**sht_info)[0]
    
    # convert alm 2D array to alm 1D array
    lmpy = len(almpy)
    alm  = libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy,lmax)
    
    return alm * 4*np.pi/(12*nside**2)


def map2alm_spin(lmax: int, mmax: int, spin: int, map_in: Array, nthreads=0, maxiter=0) -> Array:
    # expected input:
    # map_in[0,pix] = first spin component, e.g. Q
    # map_in[1,pix] = second spin component, e.g. U

    nside = int(np.sqrt(len(map_in[0]) / 12))

    hb = ducc0.healpix.Healpix_Base(nside, "RING")
    sht_info = hb.sht_info()

    almpy = ducc0.sht.adjoint_synthesis(map=map_in,spin=spin,lmax=lmax,mmax=mmax,nthreads=nthreads,**sht_info)
    almpy = np.asarray(almpy, dtype=np.complex128)

    lmpy = len(almpy[0])
    alm0 = libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy[0],lmax)
    alm1 = libcurvedsky.utils.lm_healpy2healpix(lmpy,almpy[1],lmax)

    return np.asarray([alm0, alm1]) * 4*np.pi/(12*nside**2)
