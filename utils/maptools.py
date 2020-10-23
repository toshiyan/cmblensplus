import numpy as np
import healpy as hp
import sys

import curvedsky as cs


def map_cutout( nside, lon, lat, theta, phi ):
    
    region = np.ones(12*nside**2)

    if lon[0]>lon[1]:
        sys.exit('lon is strage')
    if lat[0]>lat[1]:
        sys.exit('lat is strage')
    
    lonc = (lon[1]+lon[0])*.5
    lonw = (lon[1]-lon[0])*.5
    latc = (lat[1]+lat[0])*.5
    latw = (lat[1]-lat[0])*.5

    region[ np.abs(lonc-theta) > lonw ] = 0.
    region[ np.abs(latc-phi) > latw ] = 0.
    
    return region


def create_region( nside, lonras, latras, ascale=0. ):
    
    pix_theta, pix_phi = hp.pix2ang(nside, np.arange(12*nside**2), lonlat=True)
    region = map_cutout(nside,lonras,latras,pix_theta,pix_phi)
    
    if ascale!=0.:
        region = cs.utils.apodize(nside,region,ascale)
    
    return region


