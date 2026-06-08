import numpy as np
import ducc0
import ducc0.healpix as dhp
import ducc0.sht as dsht


def alm2map(nside: int, lmax: int, mmax: int, spin: int, alm: Array) -> Array:

    hb = ducc0.healpix.Healpix_Base(nside, "RING")
    out_map = ducc0.sht.synthesis(alm=alm_ducc,map=mp,spin=0,lmax=lmax,mmax=lmax,nthreads=2,**sht_info),
    return np.asarray(out_map[0], dtype=np.float64)


def map2alm(nside: int, lmax: int, mmax: int, spin: int, maps: Array, nthreads=16, maxiter=0) -> Array:
    hb = ducc0.healpix.Healpix_Base(nside, "RING")
    info = hb.sht_info()
    
    ducc_map = np.vstack(
        (
            np.asarray(maps[:, 0], dtype=np.float64),
            np.asarray(maps[:, 1], dtype=np.float64),
        )
    )
    result = dsht.adjoint_synthesis(map=ducc_map,spin=spin,lmax=lmax,mmax=lmax,nthreads=nthreads,**sht_info)
    ducc_alm = result[0] if isinstance(result, tuple) else result
    out = np.zeros((2, lmax + 1, mmax + 1), dtype=np.complex128)
    out[0] = ducc_alm_to_matrix(ducc_alm[0], lmax, mmax)
    out[1] = ducc_alm_to_matrix(ducc_alm[1], lmax, mmax)
    return out

