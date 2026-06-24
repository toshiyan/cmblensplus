import numpy as np
import healpy as hp
import sys

import cmblensplus.curvedsky as cs


#////////// Healpix Operation //////////#
def change_coord(m, coord):
    import healpy as hp
    """ Change coordinates of a HEALPIX map 
    taken from the following site:
      https://stackoverflow.com/questions/44443498/how-to-convert-and-save-healpy-map-to-different-coordinate-system?noredirect=1&lq=1

    Parameters
    ----------
    m : map or array of maps
      map(s) to be rotated
    coord : sequence of two character
      First character is the coordinate system of m, second character
      is the coordinate system of the output map. As in HEALPIX, allowed
      coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

    Example
    -------
    The following rotate m from galactic to equatorial coordinates.
    Notice that m can contain both temperature and polarization.
    >>>> change_coord(m, ['G', 'C'])
    """
    # Basic HEALPix parameters
    npix = m.shape[-1]
    nside = hp.npix2nside(npix)
    ang = hp.pix2ang(nside, np.arange(npix))

    # Select the coordinate transformation
    rot = hp.Rotator(coord=reversed(coord))

    # Convert the coordinates
    new_ang = rot(*ang)
    new_pix = hp.ang2pix(nside, *new_ang)

    return m[..., new_pix]


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



#////////// Rotation of Stokes Q/U map polarization angle //////////#

def qurotate(alpha,Q,U,smalla=False):
    # alpha is in unit of deg for default
    angle = alpha * np.pi/180. # deg to rad

    if smalla:
        rQ = Q - U*2*angle
        rU = Q*2*angle + U
    else:
        rQ = Q*np.cos(2*angle) - U*np.sin(2*angle)
        rU = Q*np.sin(2*angle) + U*np.cos(2*angle)
    return rQ, rU

def ebrotate(alpha,Elm,Blm,smalla=False):
    # alpha is in unit of deg for default
    
    angle = alpha * np.pi/180. # deg to rad

    if smalla:
        rElm = Elm - Blm*2*angle
        rBlm = Elm*2*angle + Blm
    else:
        rElm = Elm*np.cos(2*angle) - Blm*np.sin(2*angle)
        rBlm = Elm*np.sin(2*angle) + Blm*np.cos(2*angle)
    return rElm, rBlm




