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
    self.px1 = 0.
    self.ox1 = 0.
    self.sx1 = 0.

    self.onlydiag = False


  def x2PTE(self,simamp=1,diag=False):
    # compute chi^2 PTE of ocl using scl
    # simamp /=1 if needed
    n   = len(self.scl[:,0])
    # for real data
    mx  = np.mean(self.scl,axis=0)
    dxi = self.scl - mx
    dx0 = self.ocl - mx*simamp
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


  def x1PTE(self,simamp=1,twoside=True):
    # compute chi PTE of ocl using scl
    n   = len(self.scl[:,0])
    # for real data
    mx  = np.mean(self.scl,axis=0)
    sx  = np.std(self.scl,axis=0)
    oX1 = np.sum((self.ocl-mx*simamp)/sx)
    # for sim (exclude self rlz)
    dxi = np.array([self.scl[i,:]-np.mean(np.delete(self.scl,i,0),axis=0) for i in range(n)])
    sX1 = np.array([np.sum(dxi[i,:]/np.std(np.delete(self.scl,i,0),axis=0)) for i in range(n)])
    # output
    px1 = (sX1>oX1).sum()/np.float(n)
    if twoside:
      self.px1 = 1.-2*np.abs(px1-0.5)
    else:
      self.px1 = px1
    self.ox1 = np.around(oX1,decimals=3)
    self.sx1 = np.around(sX1,decimals=3)


  def get_amp(self,fcl=None,scale=1.,diag=False,cor=''):
    """
    Estimate the amplitude of the power spectrum

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


  def get_amp_simfix(self,ocls):
    # estimating bias in the amplitude of the power spectrum with fixed cl and covariance
    # scl is a simulated covariance and fiducial cl
    #   ocls --- mock observed cls array like [sim,bin]

    # fiducial spectrum
    fcl = np.mean(self.scl,axis=0)

    # amplitude parameters for each mock observed and simulated cl
    ampo = ocls/fcl
    amps = self.scl/fcl

    # optimal weighting evaluated from simulation
    wb, wt = opt_weight(amps,self.onlydiag)

    # amplitude estimator
    A  = np.sum(wb*ampo,axis=1)/wt
    mA = np.mean(A)
    sA = np.sqrt(np.var(A))

    return mA, sA


#////////// data analysis functions //////////#
def hist_errorbars( data, divbymax=True, xerrs=False, *args, **kwargs) :
    """Plot a histogram with error bars. Accepts any kwarg accepted by either numpy.histogram or pyplot.errorbar"""
    import matplotlib.pyplot as plt
    import inspect

    # pop off normed kwarg, since we want to handle it specially
    norm = False
    if 'normed' in kwargs.keys() :
        norm = kwargs.pop('normed')

    # retrieve the kwargs for numpy.histogram
    histkwargs = {}
    for key, value in kwargs.items():
        if key in inspect.getfullargspec(np.histogram).args :
            histkwargs[key] = value

    histvals, binedges = np.histogram( data, **histkwargs )
    a, binedges = np.histogram( data, **histkwargs)
    yerrs = np.sqrt(a)*histvals[0]/a[0]

    if norm :
        nevents = float(sum(histvals))
        binwidth = (binedges[1]-binedges[0])
        histvals = histvals/nevents/binwidth
        yerrs = yerrs/nevents/binwidth

    bincenters = (binedges[1:]+binedges[:-1])/2

    if xerrs :
        xerrs = (binedges[1]-binedges[0])/2
    else :
        xerrs = None

    # retrieve the kwargs for errorbar
    ebkwargs = {}
    for key, value in kwargs.items():
        if key in inspect.getfullargspec(plt.errorbar).args :
            histkwargs[key] = value
    if divbymax:
        histmax = np.max(histvals)
    else:
        histmax = 1.
    out = plt.errorbar(bincenters, histvals/histmax, yerrs/histmax, xerrs, fmt=".", **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')

    if 'range' in kwargs.keys() :
        plt.xlim(*kwargs['range'])

    return out


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


def get_corrcoef(scl):
    """
    Estimate correlation coefficients of the power spectrum
    """

    corr = np.corrcoef(scl,rowvar=0)
    return corr


def plot_corrcoef(b,dat,fname='',xaname='$x$',yaname='$y$',cbname='correlation coefficient',vmin=-1,vmax=1):
    import matplotlib.pyplot as plt
    x = np.linspace(b[0],b[-1],len(b))
    plt.xlim(x[0],x[-1])
    plt.ylim(x[0],x[-1])
    plt.xlabel(xaname)
    plt.ylabel(yaname)
    plt.pcolor(x,x,dat,vmin=vmin,vmax=vmax)
    cb = plt.colorbar()
    cb.set_label(cbname,labelpad=20,rotation=270)
    if fname!='': plt.savefig(fname+'.png')
    plt.show()


# Noise
def nl_white(sigma,theta,lmax,Tcmb=2.72e6):
    ac2rad = np.pi/10800.
    noise  = (sigma*ac2rad/Tcmb)**2
    L      = np.linspace(0,lmax,lmax+1)
    beam   = np.exp(L*(L+1)*(theta*ac2rad)**2/(8.*np.log(2.)))
    return noise*beam


# beam
def beam(theta,lmax):
    ac2rad = np.pi/10800.
    L      = np.linspace(0,lmax,lmax+1)
    beam   = np.exp(L*(L+1)*(theta*ac2rad)**2/(16.*np.log(2.)))
    return beam


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




