import numpy as np

# * aps amplitude histogram
class statistics:

    def __init__(self,ocl=1.,scl=1.):
        #   ocl --- observed cl array like [bin]
        #   scl --- simulated cl array like [sim,bin]
        self.ocl = ocl
        self.scl = scl
        self.A  = 0.
        self.mA = 0.
        self.sA = 0.
        self.MA = 0.
        self.p  = 0.
        self.px2 = 0.
        self.ox2 = 0.
        self.sx2 = 0.
        self.mx2 = 0.
        self.px1 = 0.
        self.ox1 = 0.
        self.sx1 = 0.
        self.mx1 = 0.
        self.onlydiag = False


    def x2PTE(self,diag=False):
        # compute chi^2 PTE of ocl using scl
        n   = len(self.scl[:,0])
        # for real data
        mx  = np.mean(self.scl,axis=0)
        dxi = self.scl - mx
        dx0 = self.ocl - mx
        cov = np.cov(dxi,rowvar=0)
        if diag: cov = np.diag(np.diag(cov))
        oX2 = np.dot(dx0,np.dot(np.linalg.inv(cov),dx0))
        # for sim (exclude self rlz)
        dxi = np.array([self.scl[i,:]-np.mean(np.delete(self.scl,i,0),axis=0) for i in range(n)])
        sX2 = np.array([np.dot(dxi[i,:],np.dot(np.linalg.inv(np.cov(np.delete(dxi,i,0),rowvar=0)),dxi[i,:])) for i in range(n)])
        # output
        self.px2 = (sX2>oX2).sum()/np.float(n)
        self.ox2 = oX2
        self.sx2 = sX2
        self.mx2 = np.mean(sX2)


    def x1PTE(self,twoside=True):
        # compute chi PTE of ocl using scl
        n   = len(self.scl[:,0])
        # for real data
        mx  = np.mean(self.scl,axis=0)
        sx  = np.std(self.scl,axis=0)
        oX1 = np.sum((self.ocl-mx)/sx)
        # for sim (exclude self rlz)
        dxi = np.array([self.scl[i,:]-np.mean(np.delete(self.scl,i,0),axis=0) for i in range(n)])
        sX1 = np.array([np.sum(dxi[i,:]/np.std(np.delete(self.scl,i,0),axis=0)) for i in range(n)])
        # output
        px1 = (sX1>oX1).sum()/np.float(n)
        if twoside:
            self.px1 = 1.-2*np.abs(px1-0.5)
        else:
            self.px1 = px1
        self.ox1 = oX1
        self.sx1 = sX1
        self.mx1 = np.mean(sX1)


    def get_amp(self,fcl=None,scale=1.,diag=False,cor='',twoside=True):
        """
        Statistics of the amplitude of the power spectrum
    
        Args:
            - statistics

        Args(optional):
            - fcl: fiducial cl used to define amplitude array like [bin]
        """

        # mean
        if fcl is None: fcl = np.mean(self.scl,axis=0)*scale
        amp = self.scl/fcl

        # covariance
        cov = np.cov(amp,rowvar=0)
        cov[np.isnan(cov)] = 0.
        if diag: cov = np.diag(np.diag(cov))

        if cor!='':
            # replace |Cor|<0.5 elements with random numbers
            Cor = np.corrcoef(amp,rowvar=0)
            c   = (np.kron(np.diag(cov),np.diag(cov).T)).reshape((len(fcl),len(fcl)))
            ran = np.random.rand(len(fcl),len(fcl))
            ran = np.tril(ran) 
            ran = ran + ran.T -np.diag(ran.diagonal())
            xy  = np.where(np.abs(Cor)<.5)
            cov[xy] += (.2*ran[xy]-.1)*c[xy]**0.5

        # coefficients
        wb = np.sum(np.linalg.inv(cov),axis=0)
        wt = np.sum(wb)

        # amplitude estimator
        A  = np.sum(wb*amp,axis=1)/wt
        mA = np.mean(A)
        sA = np.sqrt(np.var(A))

        # observed amplitude
        oA = np.sum(wb*self.ocl/fcl)/wt
        #print 'obs/sim amp', oA, mA, 'sigma amp', sA, 'ratio', oA/sA

        self.A  = A
        self.mA = mA
        self.sA = sA
        self.oA = oA
        self.p  = (A>oA).sum()/np.float(len(A))
        self.MA = np.median(A)

        if twoside:
            self.p = 1.-2*np.abs(self.p-0.5)


#////////// Statistics from MC sim //////////#
def PTEs(ocb,scb,diag=False,disp=True,x1pte=False,x2pte=True,fpt=2,comment=''):

    st = statistics(ocb,scb)

    form = '{:.'+str(fpt)+'f}'
    if comment != '':
        com = '('+comment+')'
    else:
        com = ''

    if x1pte:
        statistics.x1PTE(st)
        if disp:
            print('chi:',np.around(st.ox1,decimals=1),end=' ')
            print(', chi (sim):',np.around(st.mx1,decimals=1),end=' ')
            print(', PTE:',form.format(st.px1),com)

    if x2pte:
        statistics.x2PTE(st,diag)
        if disp:
            print('chi^2:',np.around(st.ox2,decimals=1),end=' ')
            print(', chi^2 (sim):',np.around(st.mx2,decimals=1),end=' ')
            print(', PTE:',form.format(st.px2),com)

    return st


def amplitude(ocb,scb,fcb=None,diag=False,disp=True):
    st = statistics(ocb,scb)
    statistics.get_amp(st,fcl=fcb,diag=diag)
    if disp:
        print('obs A', np.round(st.oA,3), 'mean(A)', np.round(st.mA,3), 'sigma(A)', np.round(st.sA,3), 'S/N', np.round(1./st.sA,3), 'A>oA', st.p)
    return st



def get_corrcoef(scl):
    """
    Estimate correlation coefficients of the power spectrum
    """

    corr = np.corrcoef(scl,rowvar=0)
    return corr


def get_cov(scl,fcl=None,scale=1.,diag=False,cinv=False):
    """
    Compute covariance of the power spectrum

    Args:
      - statistics

    Args(optional):
      - fcl: fiducial cl used to define amplitude array like [bin]
    """

    # mean
    if fcl is None: fcl = scale
    A = scl/fcl

    # covariance
    cov = np.cov(A,rowvar=0)
    cov[np.isnan(cov)] = 0.
    if diag: cov = np.diag(np.diag(cov))
    if cinv: cov = np.linalg.inv(cov)
    return cov


def opt_weight(x,low=-1.,diag=False):
    # optimal weighting
    #   x --- data like [sim,bin]
    cov  = np.cov(x,rowvar=0)
    if diag: cov = np.diag(np.diag(cov)) # set off-diag to zero
    cov[np.isnan(cov)] = 0.
    if low > -1.: 
        corr = get_corrcoef(x)
        cov[corr<low] = 0.
    cinv = np.linalg.inv(cov)   # inverse covariance
    wb  = np.sum(cinv,axis=0)
    wt  = np.sum(wb)           # normalization
    wi  = wb/wt # optimal weight
    return wi, wb, wt


#////////// Apodization function //////////#

def apofunc(distance,aposcale):
    # apodization window
    #   distance --- unit in rad
    #   aposcale --- scale of apodization in deg
    x = (1.-np.cos(distance))/(1.-np.cos(aposcale*np.pi/180.))
    x[x>=1.] = 1.
    y = np.sqrt(x)
    return  y - np.sin(2*np.pi*y)/(2*np.pi)


def apofunc_flat(mapsize,s,aposcale):
    # apodization window
    #   mapsize  --- unit in rad
    #   s        --- coordinates
    #   aposcale --- scale of apodization in deg
    a  = mapsize*.5
    ss = abs(s)/a
    x = (1.-ss)/(1.-aposcale)
    x[x>=1.] = 1.
    x[x<=0.] = 0.
    return  x - np.sin(2*np.pi*x)/(2*np.pi)


def window_2d(nx,ny,Dx,Dy,aposcale):
    sx = Dx/nx
    sy = Dy/ny
    xi = (np.linspace(0,nx,nx)-1.-nx*0.5)*sx
    xj = (np.linspace(0,ny,ny)-1.-ny*0.5)*sy
    Wx = apofunc_flat(Dx,xi,aposcale)
    Wy = apofunc_flat(Dy,xj,aposcale)
    return np.outer(Wx,Wy)


#////////// Likelihood function //////////#

def lnLHL(rx,fcl,icov,bi=None):
    # rx = ocb/scb
    # icov is the covariance of ocb
    bn, bn = np.shape(icov)
    if bi is None:
        gx = np.sign(rx-1.)*np.sqrt(2.*(rx-np.log(rx)-1.))
    else:
        gx = np.zeros(bn)
        gx[bi] = np.sign(rx[bi]-1.)*np.sqrt(2.*(rx[bi]-np.log(rx[bi])-1.))
    return -0.5*np.dot(gx*fcl,np.dot(icov,gx*fcl))


def lnLHLs(rx,fcl,icov,bi=None):
    # rx = ocb/scb
    # icov is the covariance of ocb
    bn, bn = np.shape(icov)
    gx = np.sign(rx-1.)*np.sqrt(2.*(rx-np.log(rx)-1.))
    return -0.5*gx*fcl[bi]*icov[bi,bi]*gx*fcl[bi]



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


#////////// Alm Operation //////////#

def ebrotate(alpha,Elm,Blm,smalla=True):
    if smalla:
        rElm = Elm - Blm*2*alpha
        rBlm = Elm*2*alpha + Blm
    else:
        rElm = Elm*np.cos(2*alpha) - Blm*np.sin(2*alpha)
        rBlm = Elm*np.sin(2*alpha) + Blm*np.cos(2*alpha)
    return rElm, rBlm



#////////// Absolute Angle Estimator //////////#

def est_absangle(oCX,sCX,oCY,sCY,fcl=1.,disp=True,diag=False,x1pte=False,x2pte=True):

    # estimate amplitude of the cross spectrum
    ocl = oCX/(oCY*2*np.pi/180.)
    scl = sCX/(sCY*2*np.pi/180.)
    st = statistics(ocl,scl)
    statistics.get_amp(st,fcl,diag=diag)

    # check PTE of observed cls
    if x1pte:
        statistics.x1PTE(st)
        print('x-PTEs of the spectrum ratio:',np.around(st.px1,decimals=3))
    if x2pte:
        statistics.x2PTE(st,diag=diag)
        print('x^2-PTEs of the spectrum ratio:',np.around(st.px2,decimals=3))

    return st


