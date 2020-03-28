
import configparser
import sys
import os
import numpy as np
import basic
import curvedsky
import misctools

if sys.version_info[:3] > (3,0):
    import pickle
elif sys.version_info[:3] > (2,5,2):
    print('use cPickle')
    import cPickle as pickle


# Define quad estimator names
class quad_fname:

    def __init__(self,pquad,qest,root,ids,cmbtag):

        # qtype is the type of mode coupling, such as lens, rot, etc
        qalm = root + pquad.qtype + '/alm/'
        qrdn = root + pquad.qtype + '/rdn0/'
        qmlm = root + pquad.qtype + '/mean/'
        qaps = root + pquad.qtype + '/aps/'

        otag = pquad.otag
        qtag = qest+'_'+cmbtag+pquad.ltag+pquad.qtagext

        # normalization and tau transfer function
        self.al   = qaps+'Al_'+qtag+'.dat'
        self.wl   = qaps+'Wl_'+qtag+'.dat'

        # N0 bias
        self.n0bs = qaps+'n0_'+qtag+'_n'+str(pquad.n0sim).zfill(3)+'.dat'

        # mean field
        self.ml   = [qmlm+'cl_'+qtag+'_'+x+'.dat' for x in ids]
        self.mfb  = [qmlm+'mlm_'+qtag+'_'+x+'.pkl' for x in ids]
        self.mfcl = qmlm+'mfcl_'+qtag+'_n'+str(pquad.mfsim).zfill(3)+'.dat'
        self.mf   = qmlm+'mfalm_'+qtag+'_n'+str(pquad.mfsim).zfill(3)+'.pkl'

        # reconstructed spectra
        self.mcls = qaps+'cl_'+qtag+'.dat'
        self.mcbs = qaps+'cl_'+qtag+otag+'.dat'
        self.ocls = qaps+'cl_'+ids[0]+'_'+qtag+'.dat'
        self.ocbs = qaps+'cl_'+ids[0]+'_'+qtag+otag+'.dat'
        self.cl   = [qaps+'rlz/cl_'+qtag+'_'+x+'.dat' for x in ids]

        # reconstructed alm and RDN0
        self.alm  = [qalm+'alm_'+qtag+'_'+x+'.pkl' for x in ids]
        self.rdn0 = [qrdn+'rdn0_'+qtag+'_n'+str(pquad.rdsim).zfill(3)+'_'+x+'.dat' for x in ids]
        self.ddn0 = [qrdn+'ddn0_'+qtag+'_'+x+'.dat' for x in ids]


class quad:

    def __init__(self,conf,qDO=None,qMV=None,qlist=None,qtype=''):

        # Define parameters for quadratic estimator computations
        self.qtype = conf.get('qtype','lens')
        if qtype!='':  self.qtype = qtype

        self.nside  = conf.getint('nside',1024)
        self.rlmin  = conf.getint('rlmin',500)
        self.rlmax  = conf.getint('rlmax',2048)
        self.oLmin  = conf.getint('oLmin',1)
        self.oLmax  = conf.getint('oLmax',2048)
        self.bn     = conf.getint('bn',30) 
        self.binspc = conf.get('binspc','')
        self.n0min  = conf.getint('n0min',1)
        self.n0max  = conf.getint('n0max',50)
        self.rdmin  = conf.getint('rdmin',1)
        self.rdmax  = conf.getint('rdmax',100)
        self.mfmin  = conf.getint('mfmin',1)
        self.mfmax  = conf.getint('mfmax',100)
        self.qtagext = conf.get('qtagext','')

        self.oL     = [self.oLmin,self.oLmax]
        self.n0sim  = self.n0max - self.n0min + 1
        self.rdsim  = self.rdmax - self.rdmin + 1
        self.mfsim  = self.mfmax - self.mfmin + 1
        self.rd4sim = conf.getboolean('rd4sim',True)  # whether RD calculation for sim

        #definition of T+P
        if qDO==None:
            self.qDO = [True,True,True,False,False,False]
        if qMV==None:
            self.qMV = ['TT','TE','EE']

        #definition of qest
        if qlist==None:
            self.qlist = ['TT','TE','EE','TB','EB','MV']
            if self.qtype=='rot':
                self.qlist = ['EB']
        else:
            self.qlist = qlist

        #cinv diag filter
        self.Fl = {}

        #multipole bins
        self.eL = np.linspace(0,self.oLmax,self.oLmax+1)
        self.bp, self.bc = basic.aps.binning(self.bn,self.oL,spc=self.binspc)

        #kappa
        self.kL = self.eL*(self.eL+1.)*.5

        #filename tags
        self.ltag = '_l'+str(self.rlmin)+'-'+str(self.rlmax)
        self.otag = '_oL'+str(self.oLmin)+'-'+str(self.oLmax)+'_b'+str(self.bn)+self.binspc

        # cmb alm mtype
        self.mtype = []
        for q in self.qlist:
            if not q[0] in self.mtype: self.mtype.append(q[0])
            if not q[1] in self.mtype: self.mtype.append(q[1])


    def fname(self,root,ids,cmbtag):
        #setup filename

        f = {}
        for q in self.qlist:
            f[q] = quad_fname(self,q,root,ids,cmbtag)

        self.f = f


    def cinvfilter(self,ocl=None):

        for m in ['T','E','B']:
            self.Fl[m] = np.zeros(self.rlmax+1)

        if ocl is None:
            for m in ['T','E','B']:
                self.Fl[m][self.rlmin:self.rlmax+1] = 1.
        else:
            for l in range(self.rlmin,self.rlmax+1):
                self.Fl['T'][l] = 1./ocl[0,l]
                self.Fl['E'][l] = 1./ocl[1,l]
                self.Fl['B'][l] = 1./ocl[2,l]


    # compute normalization
    def al(self,lcl,ocl,output=True,verbose=True,overwrite=False):
        '''
        Return normalization of the quadratic estimators
        '''

        lcl = lcl[:,:self.rlmax+1]
        ocl = ocl[:,:self.rlmax+1]

        Lmax  = self.oLmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        qlist = self.qlist
        qtype = self.qtype
        qDO   = self.qDO
        oL    = np.linspace(0,Lmax,Lmax+1)

        Ags = {}
        Acs = {}
        for q in qlist:

            if misctools.check_path(self.f[q].al,overwrite=overwrite,verbose=verbose): continue

            if qtype=='lens':
                if q=='TT': Ag, Ac = curvedsky.norm_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:],gtype='k')
                if q=='TE': Ag, Ac = curvedsky.norm_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[1,:],gtype='k')
                if q=='EE': Ag, Ac = curvedsky.norm_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],gtype='k')
                if q=='TB': Ag, Ac = curvedsky.norm_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:],gtype='k')
                if q=='EB': Ag, Ac = curvedsky.norm_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:],gtype='k')
                if q=='MV':
                    ag, ac, Wg, Wc = curvedsky.norm_lens.qall(qDO,Lmax,rlmin,rlmax,lcl,ocl,gtype='k')
                    Ag, Ac = ag[5,:], ac[5,:]

            if qtype=='rot':
                Ac = np.zeros(Lmax+1)
                if q=='TB': Ag = curvedsky.norm_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:])
                if q=='EB': Ag = curvedsky.norm_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

            if qtype=='tau':
                Ac = np.zeros(Lmax+1)
                if q=='TT': Ag = curvedsky.norm_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
                if q=='EB': Ag = curvedsky.norm_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

            # save
            if output:
                np.savetxt(self.f[q].al,np.array((oL,Ag,Ac)).T)
                if q=='MV' and qtype=='lens': 
                    for qi, qq in enumerate(['TT','TE','EE','TB','EB']): np.savetxt(self.f[qq].wl,np.array((oL,Wg[qi,:],Wc[qi,:])).T)
            else:
                Ags[q] = Ag
                Acs[q] = Ac

        if not output:
            return Ags, Acs


    def loadnorm(self):

        Ag, Ac, Wg, Wc = {}, {}, {}, {}

        # load normalization
        for q in self.qlist:
            Ag[q], Ac[q] = np.loadtxt(self.f[q].al,unpack=True,usecols=(1,2))

        # load optimal weights
        if 'MV' in self.qlist and self.qtype=='lens':
            for qi, qq in enumerate(['TT','TE','EE','TB','EB']):  Wg[qq], Wc[qq] = np.loadtxt(self.f[qq].wl,unpack=True,usecols=(1,2))

        return Ag, Ac, Wg, Wc


    def qrec(self,snmin,snmax,falm,LCl,qout=None,overwrite=False,verbose=True):
        '''
        Return quadratic estimators
        '''
        lcl = LCl[:,:self.rlmax+1]
        Lmax  = self.oLmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        qtype = self.qtype
        if qout==None:  qout = self
        rlz   = np.linspace(snmin,snmax,snmax-snmin+1,dtype=np.int)

        # load normalization and weights
        Ag, Ac, Wg, Wc = quad.loadnorm(self)

        # loop for realizations
        for i in rlz:
            
            gmv, cmv = 0., 0.

            for q in self.qlist:

                if misctools.check_path(qout.f[q].alm[i],overwrite=overwrite,verbose=verbose): continue

                if verbose:  misctools.progress(i,rlz,addtext='(qrec, '+q+')')

                if 'T' in q:  Talm = self.Fl['T'][:,None] * pickle.load(open(falm['T'][i],"rb"))[:rlmax+1,:rlmax+1]
                if 'E' in q:  Ealm = self.Fl['E'][:,None] * pickle.load(open(falm['E'][i],"rb"))[:rlmax+1,:rlmax+1]
                if 'B' in q:  Balm = self.Fl['B'][:,None] * pickle.load(open(falm['B'][i],"rb"))[:rlmax+1,:rlmax+1]

                if qtype=='lens':
                    if q=='TT':  glm, clm = curvedsky.rec_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],Talm,Talm,gtype='k',nside=nside)
                    if q=='TE':  glm, clm = curvedsky.rec_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],Talm,Ealm,gtype='k',nside=nside)
                    if q=='TB':  glm, clm = curvedsky.rec_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],Talm,Balm,gtype='k',nside=nside)
                    if q=='EE':  glm, clm = curvedsky.rec_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],Ealm,Ealm,gtype='k',nside=nside)
                    if q=='EB':  glm, clm = curvedsky.rec_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm,gtype='k',nside=nside)
                    if q=='MV':  glm, clm = gmv, cmv

                if qtype=='rot':
                    if q=='TB':  glm = curvedsky.rec_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],Talm,Balm,nside=nside)
                    if q=='EB':  glm = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm,nside=nside)
                    clm = glm*0.

                if qtype=='tau':
                    if q=='TT':  glm = curvedsky.rec_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],Talm,Talm,nside=nside)
                    if q=='EB':  glm = curvedsky.rec_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],Ealm,Balm,nside=nside)
                    clm = glm*0.

                glm *= Ag[q][:,None]
                clm *= Ac[q][:,None]
                pickle.dump((glm,clm),open(qout.f[q].alm[i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

                # T+P
                if q in self.qMV and 'MV' in self.qlist:
                    gmv += Wg[q][:,None]*glm
                    cmv += Wc[q][:,None]*clm


    def n0(self,falm,w4,LCl,overwrite=False,verbose=True):
        '''
        The N0 bias calculation
        '''

        for q in self.qlist:
            if misctools.check_path(self.f[q].n0bs,overwrite=overwrite,verbose=verbose): return

        # load normalization and weights
        Ag, Ac, Wg, Wc = quad.loadnorm(self)

        # maximum multipole of output
        lcl = LCl[:,:self.rlmax+1]
        Lmax  = self.oLmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        qlist = self.qlist
        qtype = self.qtype
        oL = np.linspace(0,Lmax,Lmax+1)
        rlz   = np.linspace(self.n0min,self.n0max,self.n0max-self.n0min+1,dtype=np.int)

        # power spectrum
        cl ={}
        for q in qlist:
            cl[q] = np.zeros((2,Lmax+1))

        # loop for realizations
        for i in rlz:

            id0 = 2*i-1
            id1 = 2*i

            gmv, cmv = 0., 0.

            for q in qlist:

                if verbose:  misctools.progress(i,rlz,addtext='(n0 bias, '+q+')')

                if q == 'MV':
                    glm, clm = gmv, cmv
                else:
                    q1, q2 = q[0], q[1]
                    if verbose:  print(q1,q2)
                    alm1 = self.Fl[q1][:,None] * pickle.load(open(falm[q1][id0],"rb"))[:rlmax+1,:rlmax+1]
                    alm2 = self.Fl[q1][:,None] * pickle.load(open(falm[q1][id1],"rb"))[:rlmax+1,:rlmax+1]

                    if q1 == q2:
                        blm1 = alm1
                        blm2 = alm2
                    else:
                        blm1 = self.Fl[q2][:,None] * pickle.load(open(falm[q2][id0],"rb"))[:rlmax+1,:rlmax+1]
                        blm2 = self.Fl[q2][:,None] * pickle.load(open(falm[q2][id1],"rb"))[:rlmax+1,:rlmax+1]
                    glm, clm = qXY(qtype,q,Lmax,rlmin,rlmax,nside,lcl,alm1,alm2,blm1,blm2)

                glm *= Ag[q][:,None]
                clm *= Ac[q][:,None]

                cl[q][0,:] += curvedsky.utils.alm2cl(Lmax,glm)/(2*w4*self.n0sim)
                cl[q][1,:] += curvedsky.utils.alm2cl(Lmax,clm)/(2*w4*self.n0sim)

                # T+P
                if q in self.qMV and 'MV' in qlist:
                    gmv += Wg[q][:,None]*glm
                    cmv += Wc[q][:,None]*clm

        for q in qlist:
            if self.n0sim > 0:
                if verbose:  print ('save N0 data')
                np.savetxt(self.f[q].n0bs,np.concatenate((oL[None,:],cl[q])).T)


    def diagrdn0(self,snmax,lcl,ocl,frcl,verbose=True):
        
        oL = np.linspace(0,self.oLmax,self.oLmax+1)
        Ag, Ac, Wg, Wc = quad.loadnorm(self)

        for i in range(snmax+1):

            if verbose:  print ('Diag-RDN0',i)

            rcl = np.loadtxt(frcl[i],unpack=True,usecols=(1,2,3,4))
            rcl[np.where(rcl==0)] = 1e30 # a large number

            # data x data
            cl = ocl**2/rcl
            Ags0, Acs0 = quad.al(self,lcl,cl,output=False)

            # (data-sim) x (data-sim)
            cl = ocl**2/(ocl-rcl)
            Ags1, Acs1 = quad.al(self,lcl,cl,output=False)

            for q in self.qlist:

                Ags0[q][np.where(Ags0[q]==0)] = 1e30
                Ags1[q][np.where(Ags1[q]==0)] = 1e30
                Acs0[q][np.where(Acs0[q]==0)] = 1e30
                Acs1[q][np.where(Acs1[q]==0)] = 1e30
                n0g = Ag[q]**2*(1./Ags0[q]-1./Ags1[q])
                n0c = Ac[q]**2*(1./Acs0[q]-1./Acs1[q])
                np.savetxt(self.f[q].drdn0[i],np.array((oL,n0g,n0c)).T)


    def rdn0(self,snmin,snmax,falm,w4,LCl,qout=None,overwrite=False,falms=None,verbose=True):
        '''
        The sim-data-mixed term of the RDN0 bias calculation
        '''

        # load normalization and weights
        Ag, Ac, Wg, Wc = quad.loadnorm(self)

        # maximum multipole of output
        lcl = LCl[:,:self.rlmax+1]
        Lmax  = self.oLmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        qlist = self.qlist
        qtype = self.qtype
        oL = np.linspace(0,Lmax,Lmax+1)
        if falms is None: falms = falm
        if qout is None:  qout = self
        rlz   = np.linspace(snmin,snmax,snmax-snmin+1,dtype=np.int)

        # load N0
        N0 = {}
        for q in qlist:
            N0[q] = np.loadtxt(self.f[q].n0bs,unpack=True,usecols=(1,2))

        # compute RDN0
        for i in rlz:

            # skip sim
            if not self.rd4sim and i!=0: 
                continue

            if verbose:  misctools.progress(i,rlz,addtext='(rdn0 bias)')

            # avoid overwriting
            Qlist = []
            for q in qlist:
                if misctools.check_path(qout.f[q].rdn0[i],overwrite=overwrite,verbose=verbose):  Qlist.append(q)
            if Qlist != []: 
                continue

            # power spectrum
            cl = {}
            for q in qlist:
                cl[q] = np.zeros((2,Lmax+1))

            # load alm
            almr = {}
            for cmb in self.mtype:
                almr[cmb] = self.Fl[cmb][:,None]*pickle.load(open(falm[cmb][i],"rb"))[:rlmax+1,:rlmax+1]


            # loop for I
            for I in range(self.rdmin,self.rdmax+1):

                gmv, cmv = 0., 0.

                # load alm
                alms = {}
                for cmb in self.mtype:
                    alms[cmb] = self.Fl[cmb][:,None]*pickle.load(open(falms[cmb][I],"rb"))[:rlmax+1,:rlmax+1]

                for q in qlist:

                    q1, q2 = q[0], q[1]

                    if I==i: continue
                    if verbose:  print(I)

                    if q=='MV':
                        glm, clm = gmv, cmv
                    else:
                        glm, clm = qXY(qtype,q,Lmax,rlmin,rlmax,nside,lcl,almr[q1],alms[q1],almr[q2],alms[q2])

                    glm *= Ag[q][:,None]
                    clm *= Ac[q][:,None]

                    cl[q][0,:] += curvedsky.utils.alm2cl(Lmax,glm)
                    cl[q][1,:] += curvedsky.utils.alm2cl(Lmax,clm)

                    # T+P
                    if q in self.qMV and 'MV' in qlist:
                        gmv += Wg[q][:,None] * glm
                        cmv += Wc[q][:,None] * clm

            if self.rdsim>0:
                
                sn = self.rdsim
                if self.rdmin<=i and i<=self.rdmax:  sn = self.rdsim-1
                
                for q in qlist:
                    cl[q] = cl[q]/(w4*sn) - N0[q]
                    if verbose:  print ('save RDN0')
                    np.savetxt(qout.f[q].rdn0[i],np.concatenate((oL[None,:],cl[q])).T)



    def mean(self,w4,overwrite=False,verbose=True):

        Lmax  = self.oLmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        qlist = self.qlist
        qtype = self.qtype
        oL = np.linspace(0,Lmax,Lmax+1)
        rlz   = np.linspace(self.mfmin,self.mfmax,self.mfmax-self.mfmin+1,dtype=np.int)

        for q in qlist:

            if misctools.check_path(self.f[q].mf,overwrite=overwrite,verbose=verbose): continue

            if verbose:  print('compute mean field')
            mfg = 0.
            mfc = 0.
            for I in rlz:
                if verbose:  misctools.progress(I,rlz,addtext='(mean)')
                mfgi, mfci = pickle.load(open(self.f[q].alm[I],"rb"))
                mfg += mfgi/self.mfsim
                mfc += mfci/self.mfsim

            if verbose:  print('save to file')
            pickle.dump((mfg,mfc),open(self.f[q].mf,"wb"),protocol=pickle.HIGHEST_PROTOCOL)

            # compute mf cls
            if verbose:  print('cl for mean field bias')
            cl = np.zeros((2,Lmax+1))
            cl[0,:] = curvedsky.utils.alm2cl(Lmax,mfg)/w4
            cl[1,:] = curvedsky.utils.alm2cl(Lmax,mfc)/w4
            np.savetxt(self.f[q].mfcl,np.concatenate((oL[None,:],cl)).T)



    def mean_rlz(self,snmin,snmax,w4,overwrite=False,verbose=True):

        Lmax  = self.oLmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        qlist = self.qlist
        qtype = self.qtype
        oL = np.linspace(0,Lmax,Lmax+1)
        rlz   = np.linspace(snmin,snmax,snmax-snmin+1,dtype=np.int)

        for q in qlist:

            glm = np.zeros((self.snmf,Lmax+1,Lmax+1),dtype=np.complex)
            clm = np.zeros((self.snmf,Lmax+1,Lmax+1),dtype=np.complex)

            for I in range(1,self.snmf+1):
                if verbose:  misctools.progress(I,np.linspace(1,self.snmf,self.snmf),addtext='(load reconstructed alms, '+q+')')
                glm[I-1,:,:], clm[I-1,:,:] = pickle.load(open(self.f[q].alm[I],"rb"))

            for i in rlz:

                if misctools.check_path(self.f[q].mfb[i],overwrite=overwrite,verbose=verbose): continue

                if verbose:  misctools.progress(i,rlz,addtext='(compute mean field for each rlz)')
                mfg = np.average(glm,axis=0)
                mfc = np.average(clm,axis=0)
                if i!=0: 
                    mfg -= glm[i-1,:,:]/self.snmf
                    mfc -= clm[i-1,:,:]/self.snmf
                    mfg *= self.snmf/(self.snmf-1.)
                    mfc *= self.snmf/(self.snmf-1.)

                if verbose:  print('save to file')
                pickle.dump((mfg,mfc),open(self.f[q].mfb[i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)


            # compute mf cls
            if verbose:  print('cl for mean field bias')

            for i in rlz:

                if misctools.check_path(self.f[q].ml[i],overwrite=overwrite,verbose=verbose): continue

                mfg, mfc = pickle.load(open(self.f[q].mfb[i],"rb"))

                cl = np.zeros((2,Lmax+1))
                cl[0,:] = curvedsky.utils.alm2cl(Lmax,mfg) / w4
                cl[1,:] = curvedsky.utils.alm2cl(Lmax,mfc) / w4
                np.savetxt(self.f[q].ml[i],np.concatenate((oL[None,:],cl)).T)



def qXY(qtype,qcomb,Lmax,rlmin,rlmax,nside,lcl,alm1,alm2,blm1,blm2):
        # for N0 and RDN0 estimates

        if qtype=='lens':
            if qcomb=='TT':
                glm1, clm1 = curvedsky.rec_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],alm1,alm2,gtype='k',nside=nside)
                glm2, clm2 = curvedsky.rec_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],alm2,alm1,gtype='k',nside=nside)
            if qcomb=='TE':
                glm1, clm1 = curvedsky.rec_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],alm1,blm2,gtype='k',nside=nside)
                glm2, clm2 = curvedsky.rec_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],alm2,blm1,gtype='k',nside=nside)
            if qcomb=='TB':
                glm1, clm1 = curvedsky.rec_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm1,blm2,gtype='k',nside=nside)
                glm2, clm2 = curvedsky.rec_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm2,blm1,gtype='k',nside=nside)
            if qcomb=='EE':
                glm1, clm1 = curvedsky.rec_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],alm1,alm2,gtype='k',nside=nside)
                glm2, clm2 = curvedsky.rec_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],alm2,alm1,gtype='k',nside=nside)
            if qcomb=='EB':
                glm1, clm1 = curvedsky.rec_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm1,blm2,gtype='k',nside=nside)
                glm2, clm2 = curvedsky.rec_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm2,blm1,gtype='k',nside=nside)

            return glm1+glm2, clm1+clm2

        if qtype=='rot':
            if qcomb=='TB':
                rlm1 = curvedsky.rec_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm1,blm2,nside=nside)
                rlm2 = curvedsky.rec_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm2,blm1,nside=nside)
            if qcomb=='EB':
                rlm1 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm1,blm2,nside=nside)
                rlm2 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm2,blm1,nside=nside)

            return rlm1+rlm2, (rlm1+rlm2)*0.

        if qtype=='tau':
            if qcomb=='TT':
                rlm1 = curvedsky.rec_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],alm1,blm2,nside=nside)
                rlm2 = curvedsky.rec_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],alm2,blm1,nside=nside)
            if qcomb=='EB':
                rlm1 = curvedsky.rec_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm1,blm2,nside=nside)
                rlm2 = curvedsky.rec_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm2,blm1,nside=nside)

            return rlm1+rlm2, (rlm1+rlm2)*0.



def n0x(qx,qd0,qd1,falm,fblm,w4,lcl):
    '''
    N0 for a^1 x a^2
    '''

    # load normalization and weights
    Ag0, Ac0, Wg0, Wc0 = quad.loadnorm(qd0)
    Ag1, Ac1, Wg1, Wc1 = quad.loadnorm(qd1)

    # maximum multipole of output
    lcl = lcl[:,:qx.rlmax+1]
    Lmax  = qx.oLmax
    rlmin = qx.rlmin
    rlmax = qx.rlmax
    nside = qx.nside
    qlist = qx.qlist
    qtype = qx.qtype
    oL = np.linspace(0,Lmax,Lmax+1)

    # power spectrum
    cl ={}
    for q in qlist:
        cl[q] = np.zeros((2,Lmax+1))

    # loop for realizations
    for i in range(qx.snn0):
        print (2*i+1, 2*i+2)

        gmv, cmv = 0., 0.

        for q in qlist:

            if q == 'MV':
                glm, clm = gmv, cmv
            else:
                q1, q2 = q[0], q[1]
                print(q1,q2)
                alm1 = qd0.Fl[q1][:,None] * pickle.load(open(falm[q1][2*i+1],"rb"))[:rlmax+1,:rlmax+1]
                alm2 = qd1.Fl[q1][:,None] * pickle.load(open(fblm[q1][2*i+2],"rb"))[:rlmax+1,:rlmax+1]
                if q1 == q2:
                    blm1 = alm1
                    blm2 = alm2
                else:
                    blm1 = qd0.Fl[q2][:,None] * pickle.load(open(falm[q2][2*i+1],"rb"))[:rlmax+1,:rlmax+1]
                    blm2 = qd1.Fl[q2][:,None] * pickle.load(open(fblm[q2][2*i+2],"rb"))[:rlmax+1,:rlmax+1]
                glm, clm = qXY(qtype,q,Lmax,rlmin,rlmax,nside,lcl,alm1,alm2,blm1,blm2)

            cl[q][0,:] += Ag0[q]*Ag1[q]*curvedsky.utils.alm2cl(Lmax,glm)/(2*w4*qx.snn0)
            cl[q][1,:] += Ac0[q]*Ac1[q]*curvedsky.utils.alm2cl(Lmax,clm)/(2*w4*qx.snn0)

            # T+P
            if q in qx.qMV and 'MV' in qlist:
                gmv += Wg[q][:,None]*glm
                cmv += Wc[q][:,None]*clm

    for q in qlist:
        if qx.snn0 > 0:
            print ('save N0 data')
            np.savetxt(qx.f[q].n0bl,np.concatenate((oL[None,:],cl[q])).T)



def rdn0x(qx,qd0,qd1,snmin,snmax,falm,fblm,w4,lcl):
    '''
    The sim-data-mixed term of the RDN0 bias calculation
    '''

    Ag0, Ac0, Wg0, Wc0 = quad.loadnorm(qd0)
    Ag1, Ac1, Wg1, Wc1 = quad.loadnorm(qd1)

    # maximum multipole of output
    lcl = lcl[:,:qx.rlmax+1]
    Lmax  = qx.oLmax
    rlmin = qx.rlmin
    rlmax = qx.rlmax
    nside = qx.nside
    qlist = qx.qlist
    qtype = qx.qtype
    oL = np.linspace(0,Lmax,Lmax+1)

    # load N0
    N0 = {}
    for q in qlist:
        N0[q] = np.loadtxt(qx.f[q].n0bl,unpack=True,usecols=(1,2))

    # compute RDN0
    for i in range(snmin,snmax+1):
        print(i)

        # power spectrum
        cl = {}
        for q in qlist:
            cl[q] = np.zeros((1,Lmax+1))

        # load alm
        almr = {}
        blmr = {}
        for cmb in qx.mtype:
            almr[cmb] = qd0.Fl[cmb][:,None]*pickle.load(open(falm[cmb][i],"rb"))[:rlmax+1,:rlmax+1]
            blmr[cmb] = qd1.Fl[cmb][:,None]*pickle.load(open(fblm[cmb][i],"rb"))[:rlmax+1,:rlmax+1]

        # loop for I
        for I in range(1,qx.snrd+1):

            gmv, cmv = 0., 0.

            # load alm
            alms = {}
            blms = {}
            for cmb in qx.mtype:
                alms[cmb] = qd0.Fl[cmb][:,None]*pickle.load(open(falm[cmb][I],"rb"))[:rlmax+1,:rlmax+1]
                blms[cmb] = qd1.Fl[cmb][:,None]*pickle.load(open(fblm[cmb][I],"rb"))[:rlmax+1,:rlmax+1]

            for q in qlist:

                q1, q2 = q[0], q[1]

                if I==i: continue
                print(I)

                if q=='MV':
                    glm, clm = gmv, cmv
                else:
                    glm0, glm1 = qXYx(qtype,q,Lmax,rlmin,rlmax,nside,lcl,almr,blmr,alms,blms)

                cl[q][0,:] += Ag0[q]*Ag1[q]*curvedsky.utils.alm2cl(Lmax,glm0,glm1)
                #cl[q][1,:] += Ac0[q]*Ac1[q]*curvedsky.utils.alm2cl(Lmax,clm)

                # T+P
                if q in qx.qMV and 'MV' in qlist:
                    gmv += (Wg0[q][:,None]*Wg1[q][:,None])**0.5 * glm
                    #cmv += (Wc0[q][:,None]*Wc1[q][:,None])**0.5 * clm

        if qx.snrd>0:
            if i==0:  sn = qx.snrd
            if i!=0:  sn = qx.snrd-1
            for q in qlist:
                cl[q] = cl[q]/(w4*sn) - N0[q]
                print ('save RDN0')
                np.savetxt(qx.f[q].rdn0[i],np.concatenate((oL[None,:],cl[q])).T)


def qXYx(qtype,qcomb,Lmax,rlmin,rlmax,nside,lcl,almr,blmr,alms,blms):

    q1, q2 = qcomb[0], qcomb[1]

    if qtype=='rot':
        if qcomb=='EB':
            a01 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],almr[q1],alms[q2],nside=nside)
            a10 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],alms[q1],almr[q2],nside=nside)
            b01 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],blmr[q1],blms[q2],nside=nside)
            b10 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],blms[q1],blmr[q2],nside=nside)

        return a01+a10, b01+b10

