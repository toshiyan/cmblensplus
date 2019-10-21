# Flatsky lensing normalization from an assymmetric quadratic estimator
import numpy as np
from orphics import stats
from pixell import enmap

def getspec(f, lmin = 50, lmax = 4000, deltal = 20):
    p2d = enmap.read_map(f)
    shape,wcs = p2d.shape, p2d.wcs
    bin_edges = np.arange(lmin, lmax, deltal)
    modlmap = enmap.modlmap(shape, wcs)
    binner = stats.bin2D(modlmap, bin_edges)
    cents, p1d = binner.bin(p2d)
    return cents, p1d

# load unlensed and lensed Cls
bc, n0 = getspec('tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint_noise.fits')
bc, nx = getspec('tilec_single_tile_deep56_cmb_deprojects_comptony_map_v1.0.0_rc_joint_cross_noise.fits')
bc, n1 = getspec('tilec_single_tile_deep56_cmb_deprojects_comptony_map_v1.0.0_rc_joint_noise.fits')
np.savetxt('nl.dat',np.array((bc,n0,n1,nx)).T,fmt='%8.6e')


