import numpy as np
import basic

#//// Delensing ////#
def resbb_multitracer(olmax,elmin,elmax,klmin,klmax,vec,cov,kk,lEE,pp,WE):
    
    Lmax = len(vec[0,:]) - 1

    ava = np.zeros(Lmax+1)
    for L in range(Lmax+1):
        if np.sum(vec[:,L]) == 0.: continue
        ava[L] = np.dot(vec[:,L],np.dot(np.linalg.inv(cov[:,:,L]),vec[:,L])) / kk[L]
    
    tbb = basic.delens.lintemplate(olmax,elmin,elmax,klmin,klmax,lEE[:elmax+1],pp[:klmax+1],WE[:elmax+1],ava[:klmax+1])
    return tbb

    
def multitracer_weights(cov, vec, Lmin=2):
    '''
    Calculate Eq.(A8) of http://arxiv.org/abs/1502.05356 or Eq.(12) of https://arxiv.org/pdf/1511.04653.pdf
    '''
    n    = len(cov[:,0,0]) # assuming cov[n,n,Lmax]
    Lmax = len(cov[0,0,:]) - 1

    icov = np.zeros(cov.shape)
    for L in range(Lmax+1):
        try:
            icov[:,:,L] = np.linalg.inv(cov[:,:,L])
        except:
            pass

    weight = np.zeros((n,Lmax+1))
    for index in range(n):
        for L in range(Lmin,Lmax+1):
            weight[index,L] = np.dot(vec[:,l],cov[index,:,L])

    return weight


