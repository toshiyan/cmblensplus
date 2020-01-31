import numpy as np
import matplotlib.pyplot as plt


def plot_corr(dat,bp,zmin=-1,zmax=1,spc='',fname='',xlab='',ylab='',clab='corr. coeff.'):
    """
    Plot correlation coefficient of the input data [rlz,x]
    """
    #bn  = np.shape(dat)[1]
    cor = np.corrcoef(dat,rowvar=0)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if spc == '':
        x = bp
    if spc == 'p2':
        x = np.sqrt(bp)
        xs = np.array([20,50,100,200,500,1000,2000])
        plt.xticks(np.sqrt(xs),xs)

    plt.pcolor(x,x,cor,vmin=zmin,vmax=zmax)
    cb = plt.colorbar()
    cb.set_label(clab,labelpad=20,rotation=270)
    if fname!='': plt.savefig(fname+'.png')
    plt.show()


def plot_1dstyle(spc='',frac=False,xlab='$L$',ylab='$C_L$',xmin=20,xmax=2048,ymin=None,ymax=None,xlog=False,ylog=False,xylog=False,yticks=None,grid=False,fsize=None,xlabsize=14,ylabsize=14,xlabloc=None,ylabloc=None):
    """
    Start to define plot environment for 1D function
    frac --- for fractional difference (add y=0,-1,1 lines)
    """

    params = {
        'axes.labelsize': 8,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'text.usetex': False,
        }
    plt.rcParams.update(params)

    if fsize is not None:
        plt.rcParams["figure.figsize"] = (fsize[0],fsize[1])

    plt.grid(grid)
    
    if xlog or xylog: plt.xscale('log')
    if ylog or xylog: plt.yscale('log')
    
    if frac:
        if ymin is None: ymin = -5.
        if ymax is None: ymax = 5.
        plt.axhline(0,ls='--',color='k')
        plt.axhline(1,ls='--',color='k',lw=.5)
        plt.axhline(-1,ls='--',color='k',lw=.5)
    
    if yticks is not None:
        plt.yticks(yticks)

    if xlabloc is not None:
        ax = plt.gca()
        ax.xaxis.set_label_coords(xlabloc[0],xlabloc[1])

    if ylabloc is not None:
        ax = plt.gca()
        ax.yaxis.set_label_coords(ylabloc[0],ylabloc[1])

    plt.xlabel(xlab,fontsize=xlabsize)
    plt.ylabel(ylab,fontsize=ylabsize)

    if ymin is not None and ymax is not None:
        plt.ylim(ymin,ymax)
    
    if spc == '':
        plt.xlim(xmin,xmax)

    if spc == 'p2':
        plt.xlim(np.sqrt(xmin),np.sqrt(xmax))
        xs = np.array([20,50,100,200,500,1000,2000])
        plt.xticks(np.sqrt(xs),xs)
    

def hist_errorbars( data, ymin=None, ymax=None, divbymax=True, xerrs=False, *args, **kwargs) :
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
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = 1
    plt.ylim(ymin,ymax)
    out = plt.errorbar(bincenters, histvals/histmax, yerrs/histmax, xerrs, fmt=".", **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')

    if 'range' in kwargs.keys() :
        plt.xlim(*kwargs['range'])

    return out





