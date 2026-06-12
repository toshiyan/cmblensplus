
import configparser
import sys
import os
import numpy as np
import pickle
import tqdm
import multiprocessing

# cmblensplus
import cmblensplus.curvedsky as curvedsky

# local
from . import misctools


def set_mtype(qlist):  
    # set CMB alm mtype (T, E, B) from quadratic estimator combinations

    mtype = []
    for q in qlist:
        if q not in ['TT','TE','EE','TB','EB','BB','MV']:
            sys.exit('invalid quadratic combination is specified')
        if q == 'MV': continue
        if not q[0] in mtype: mtype.append(q[0])
        if not q[1] in mtype: mtype.append(q[1])

    return mtype                



# Define quadratic estimator names
class quad_fname():

    def __init__(self,qobj,qest):

        # qtype is the type of mode coupling, such as lens, rot, etc
        qalm = qobj.qalm 
        qrdn = qobj.qrdn
        qmlm = qobj.qmlm
        qaps = qobj.qaps

        qtag = qest + '_' + qobj.cmbtag + qobj.bhe_tag + qobj.ltag + qobj.qtagext

        # normalization and tau transfer function
        self.al   = qaps+'Al_'+qtag+'.dat'
        self.wl   = qaps+'Wl_'+qtag+'.dat'

        # N0 bias
        self.n0bs = qaps+'n0_'+qtag+'_n'+str(qobj.n0sim).zfill(3)+'.dat'

        # N0 bias for each rlz pair
        self.n0bs_rlz = [ qrdn+'n0_'+qtag+'_'+x+'.dat' for x in qobj.ids ]
        
        # mean field
        self.mfcl  = [qmlm+'cl_'+qtag+'_'+x+'.dat' for x in qobj.ids]
        self.mfalm = [qmlm+'mlm_'+qtag+'_'+x+'.pkl' for x in qobj.ids]
        self.MFcl  = qmlm+'mfcl_'+qtag+'_n'+str(qobj.mfsim).zfill(3)+'.dat'
        self.MFalm = qmlm+'mfalm_'+qtag+'_n'+str(qobj.mfsim).zfill(3)+'.pkl'

        # reconstructed spectra
        self.mcls = qaps+'cl_'+qtag+'.dat'
        self.ocls = qaps+'cl_'+qobj.ids[0]+'_'+qtag+'.dat'
        self.cl   = [qobj.qcls+'cl_'+qtag+'_'+x+'.dat' for x in qobj.ids]

        # reconstructed alm and RDN0
        self.alm  = [qalm+'alm_'+qtag+'_'+x+'.pkl' for x in qobj.ids]
        self.rdn0 = [qrdn+'rdn0_'+qtag+'_n'+str(qobj.rdsim).zfill(3)+'_'+x+'.dat' for x in qobj.ids]

        # RDN0 for each sim rlz
        self.rdn0_rlz = [ {I: qrdn+'rdn0_'+qtag+'_rd'+str(I).zfill(3)+'_'+x+'.dat' for I in qobj.rdids } for x in qobj.ids ]

        # diagonal RDN0
        self.drdn0 = [qrdn+'drdn0_'+qtag+'_'+x+'.dat' for x in qobj.ids]

        # reconstruction noise variance map
        self.nkmap = qalm+'nkmap_'+qtag+'.pkl'

        # additional alms
        self.walm = [qalm+'walm_'+qtag+'_'+x+'.pkl' for x in qobj.ids]
        


'''
def cmb_data(qobj,lcl=None,ocl=None,ifl=None,falm='',stag=''):

    # determined by CMB data to be used for reconstruction

    # Cl for filter, obs Cl and alm files
    qobj.ifl = ifl
    if lcl is not None:  qobj.lcl = lcl[:,:qobj.rlmax+1]
    if ocl is not None:  qobj.ocl = ocl[:,:qobj.rlmax+1]
    if ifl is not None:  qobj.ifl = ifl[:,:qobj.rlmax+1]
    qobj.falm = falm

    #cinv diag filter
    qobj.Fl = { m: np.zeros(qobj.rlmax+1) for m in ['T','E','B'] }
        
    # tag for cmb data given by hand
    qobj.cmbtag = stag
'''


def setup_bhe(qobj,bhe): # bhe types
    
    if bhe is None:
        qobj.bhe_do   = False
        qobj.bhe_list = []
    else:
        qobj.bhe_do   = True
        qobj.bhe_list = bhe
        if not isinstance(bhe,list):
            sys.exit('bhe should be a list')
        
    if qobj.qtype in qobj.bhe_list:
        sys.exit('qtype is included in bhe to be deprojected, please remove')

    if bhe is None:
        qobj.bhe_tag = ''
    else:
        qobj.bhe_tag = '_'+'-'.join(['bh']+qobj.bhe_list)

        
class quad(): 

    def __init__(self, 
                 rlz=None, ids=[], 
                 lcl=None, ocl=None, ifl=None, 
                 falm='',
                 olmax=2048, rlmin=500, rlmax=3000, nside=2048,
                 n0min=1, n0max=50, output_n0bs_rlz=False,
                 rdmin=1, rdmax=100, rd4sim=False, 
                 mfmin=1, mfmax=100,
                 qDO=None, qMV=None, qlist=None, qtype='', wn=None, bhe=None, 
                 overwrite=False, verbose=True,
                 stag='', root='', qtagext=''
                ):
        
        '''
        lcl[0:4,0:rlmax] : Theory CMB cls used for the estimator and normalization
        ocl[0:3,0:rlmax] : A best estimate of observd CMB cls (sum of signal and noise spectra) for TT(0), EE(1), BB(2). 
                           This will be used for the normalization and diagonal inverse variance filtering if necessary
        ifl[0:3,0:rlmax] : Diagonal inverse variance filter (0: 1/TT, 1: 1/EE, 2: 1/BB). If ifl=None, the code assumes that the input alms are already filtered. 
        falm[m][i]       : Filenames for m = T, E, or B and i = ith realization in rlz. For example, falm = {m: ['test'+str(i)+'.dat' for i in range(100)] for m in ['T','E','B']}
        '''

        #//// get parameters ////#
        conf = misctools.load_config('QUADREC')

        # define parameters for quadratic estimator computations
        self.qtype = conf.get('qtype','lens')
        if qtype!='':  self.qtype = qtype

        self.nside  = conf.getint('nside',nside)
        self.rlmin  = conf.getint('rlmin',rlmin)
        self.rlmax  = conf.getint('rlmax',rlmax)

        if 3*self.nside < self.rlmax:  print('Warning: nside^rec would be too small for reconstruction')

        self.olmin  = conf.getint('olmin',1)
        self.olmax  = conf.getint('olmax',olmax)

        # iteration
        if rlz is None:
            self.rlz = range(len(ids))
        else:
            self.rlz = rlz
        self.ids = ids
        
        # start, stop rlz of N0 bias
        self.n0min  = conf.getint('n0min',n0min)
        self.n0max  = conf.getint('n0max',n0max)
        self.n0sim  = self.n0max - self.n0min + 1
        self.n0rlz  = np.linspace(self.n0min,self.n0max,self.n0sim,dtype=np.int32)
        self.output_n0bs_rlz = output_n0bs_rlz

        # start, stop rlz of RDN0 bias
        self.rd4sim = conf.getboolean('rd4sim',rd4sim)  # whether RDN0 calculation for sim
        self.rdmin  = conf.getint('rdmin',rdmin) # sim index minimum for RDN0
        self.rdmax  = conf.getint('rdmax',rdmax) # sim index maximum for RDN0
        self.rdids  = np.arange(self.rdmin,self.rdmax+1)  # index for RDN0 sim
        self.rdsim  = self.rdmax - self.rdmin + 1 # number of RDN0 sim
        
        # start, stop rlz of mean-field
        self.mfmin  = conf.getint('mfmin',mfmin)
        self.mfmax  = conf.getint('mfmax',mfmax)
        self.mfsim  = self.mfmax - self.mfmin + 1

        # external tag
        self.qtagext = conf.get('qtagext',qtagext)

        self.oL     = [self.olmin,self.olmax] # minimum/maximum output multipoles

        # rlz
        self.mfrlz = np.linspace(self.mfmin,self.mfmax,self.mfsim,dtype=np.int32)

        # Input CMB stuff: Cl for filter, obs Cl and alm data files
        #cmb_data(self,lcl=lcl,ocl=ocl,ifl=ifl,falm=falm,stag=stag)
        self.ifl = ifl
        if lcl is not None:  self.lcl  = lcl[:,:self.rlmax+1]
        if ocl is not None:  self.ocl  = ocl[:,:self.rlmax+1]
        if ifl is not None:  self.ifl  = ifl[:,:self.rlmax+1]
        self.falm = falm

        #cinv diag filter
        self.Fl = { m: np.zeros(self.rlmax+1) for m in ['T','E','B'] }
   
        # tag for cmb data given by hand
        self.cmbtag = stag

        # definition of T+P
        if qDO is None:
            self.qDO = [True,True,True,False,False,False]
        else:
            self.qDO = qDO
            
        # estimators for MV
        if qMV is None:
            self.qMV = ['TT','TE','EE']
        else:
            self.qMV = qMV

        #definition of qlist
        if qlist is None:
            self.qlist = ['TT','TE','EE','TB','EB','MV']
            if self.qtype=='rot':
                self.qlist = ['EB']
        else:
            self.qlist = qlist.copy()

        # window normalization correction
        if wn is None:
            self.wn = np.ones(5)
        else:
            self.wn = wn

        #multipole bins
        self.l = np.linspace(0,self.olmax,self.olmax+1)

        #kappa
        self.kL = self.l*(self.l+1.)*.5

        #filename tags
        self.ltag = '_l'+str(self.rlmin)+'-'+str(self.rlmax)

        # cmb alm mtype
        self.mtype = set_mtype(self.qlist)
        
        #//// Bias Herdened Estimators ////
        setup_bhe(self,bhe)

        #//// File and Directory ////#
        self.root = root
        self.qalm = self.root + self.qtype + '/alm/'     # reconstructed alms
        self.qrdn = self.root + self.qtype + '/rdn0/'    # RDN0
        self.qmlm = self.root + self.qtype + '/mean/'    # Mean field
        self.qaps = self.root + self.qtype + '/aps/'     # Cls of reconstructed alm
        self.qcls = self.root + self.qtype + '/aps/rlz/' # Cls for each realization

        #setup filename
        self.f = {q: quad_fname(self,q) for q in ['TT','TE','TB','EE','EB','BB','MV']}

        #//// Misc ////#
        self.overwrite = overwrite
        self.verbose = verbose


    def create_directory(self):
        # create directories if they do not exist
        for directory in [self.qalm,self.qrdn,self.qmlm,self.qaps,self.qcls]:
            misctools.create_directory(directory,verbose=False)


    def cinvfilter(self,mids={'T':0,'E':1,'B':2}):

        if self.ifl is None:
            for m in ['T','E','B']:
                self.Fl[m][self.rlmin:self.rlmax+1] = 1.
        else:
            ocl = self.ifl.copy()
            for l in range(self.rlmin,self.rlmax+1):
                for m in self.mtype:
                    i = mids[m]
                    self.Fl[m][l] = 1./ocl[i,l]
                    if ocl[i,l] <=0.:  
                        sys.exit(m+' inverse-filter is zero: ocl='+str(ocl[i,l])+' at l='+str(l))


    # compute normalization
    def al(self,ocls=None,output=True,store=True,val_return=False,gtype='k'):
        '''
        Return normalization of the quadratic estimators
        '''

        ocl = self.ocl
        lcl = self.lcl
        if ocls is not None:  ocl = ocls

        Lmax  = self.olmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        #oL    = np.linspace(0,Lmax,Lmax+1)
        
        bhe_c = {}

        if store:
            self.Ag, self.Ac, self.Wg, self.Wc = {}, {}, {}, {}
            self.bhe_c = {}
        
        if val_return:
            Ags, Acs = {}, {}
        
        for q in self.qlist:

            if output and misctools.check_path(self.f[q].al,overwrite=self.overwrite,verbose=self.verbose): continue

            Al, At, As = 0., 0., 0. # for BHE
            
            if self.qtype=='lens' or 'lens' in self.bhe_list:
                if q=='TT': Ag, Ac = curvedsky.norm_quad.qtt('lens',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:],lfac=gtype)
                if q=='TE': Ag, Ac = curvedsky.norm_quad.qte('lens',Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[1,:],lfac=gtype)
                if q=='EE': Ag, Ac = curvedsky.norm_quad.qee('lens',Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],lfac=gtype)
                if q=='TB': Ag, Ac = curvedsky.norm_quad.qtb('lens',Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:],lfac=gtype)
                if q=='EB': Ag, Ac = curvedsky.norm_quad.qeb('lens',Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:],lfac=gtype)
                if q=='MV':
                    ag, ac, Wg, Wc = curvedsky.norm_quad.qall('lens',self.qDO,Lmax,rlmin,rlmax,lcl,ocl,lfac=gtype)
                    Ag, Ac = ag[5,:], ac[5,:]
                Al = Ag.copy() # for BHE

            if self.qtype=='rot':
                if q=='TB': Ag, Ac = curvedsky.norm_quad.qtb('rot',Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:])
                if q=='EB': Ag, Ac = curvedsky.norm_quad.qeb('rot',Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])

            if self.qtype=='tau' or 'tau' in self.bhe_list:
                if q=='TT': Ag, Ac = curvedsky.norm_quad.qtt('amp',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
                if q=='EB': Ag, Ac = curvedsky.norm_quad.qeb('amp',Lmax,rlmin,rlmax,lcl[1,:],ocl[1,:],ocl[2,:])
                At = Ag.copy() # for BHE

            if self.qtype=='src' or 'src' in self.bhe_list:
                if q=='TT': Ag, Ac = curvedsky.norm_quad.qtt('src',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])
                As = Ag.copy() # for BHE
                    
            if self.qtype=='ilens' or 'ilens' in self.bhe_list:
                if q=='TE': Ag, Ac = curvedsky.norm_imag.qte('lens',Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[1,:],lfac=gtype)
                if q=='EE': Ag, Ac = curvedsky.norm_imag.qee('lens',Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],ocl[1,:],lfac=gtype)
                if q=='TB': Ag, Ac = curvedsky.norm_imag.qtb('lens',Lmax,rlmin,rlmax,lcl[3,:],ocl[0,:],ocl[2,:],lfac=gtype)
                if q=='EB': Ag, Ac = curvedsky.norm_imag.qeb('lens',Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],ocl[1,:],ocl[2,:],lfac=gtype)
                if q=='BB': Ag, Ac = curvedsky.norm_imag.qbb('lens',Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],ocl[2,:],lfac=gtype)
                if q=='MV':
                    ag, ac, Wg, Wc = curvedsky.norm_imag.qall('lens',self.qDO,Lmax,rlmin,rlmax,lcl,ocl,lfac=gtype)
                    Ag, Ac = ag[5,:], ac[5,:]
                Al = Ag.copy() # for BHE

            #//// Bias-hardened estimator (cross response) ////#
            # Currently, only TT is supported
            Rlt, Rls, Rts = 0., 0., 0.
            
            if self.qtype == 'lens':
                if q == 'TT':
                    if 'tau' in self.bhe_list:
                        Rlt = curvedsky.norm_quad.xtt('lensamp',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:],lfac=gtype)
                    if 'src' in self.bhe_list:
                        Rls = curvedsky.norm_quad.xtt('lenssrc',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:],lfac=gtype)
                    if 'src' in self.bhe_list and 'tau' in self.bhe_list:  
                        Rts = curvedsky.norm_quad.xtt('ampsrc',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])

            if self.qtype == 'tau':
                if q == 'TT':
                    if 'lens' in self.bhe_list:
                        Rlt = curvedsky.norm_quad.xtt('lensamp',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:],lfac=gtype)
                    if 'src' in self.bhe_list and 'lens' in self.bhe_list:  
                        Rls = curvedsky.norm_quad.xtt('lenssrc',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:],lfac=gtype)
                    if 'src' in self.bhe_list:
                        Rts = curvedsky.norm_quad.xtt('ampsrc',Lmax,rlmin,rlmax,lcl[0,:],ocl[0,:])

            # Denominator
            DetR = 1 - Al*As*Rls**2 - Al*At*Rlt**2 - At*As*Rts**2 + 2.*Al*At*As*Rlt*Rls*Rts

            # Corrected normalization (to be multiplied to the unnormalized estimators)
            bhe_c[q] = {}

            if self.qtype == 'lens' and self.bhe_do:
                bhe_c[q]['lens']  = ( 1. - At*As*Rts**2 ) / DetR * Al
                bhe_c[q]['tau'] = ( Rls*As*Rts - Rlt ) / DetR * At * Al
                bhe_c[q]['src']  = ( Rlt*At*Rts - Rls ) / DetR * As * Al
                Ag = bhe_c[q]['lens']

            if self.qtype == 'tau' and self.bhe_do:
                bhe_c[q]['tau']  = ( 1. - Al*As*Rls**2 ) / DetR * At
                bhe_c[q]['lens'] = ( Rts*As*Rls - Rlt ) / DetR * Al * At
                bhe_c[q]['src']  = ( Rlt*Al*Rls - Rts ) / DetR * As * At
                Ag = bhe_c[q]['tau']
                Ac = Ag*.0

            #//// save ////#
            # output to disk
            if output:
                np.savetxt(self.f[q].al,np.array((self.l,Ag,Ac)).T)
                if q=='MV' and self.qtype=='lens':
                    for qi, qq in enumerate(self.qMV):  
                        np.savetxt(self.f[qq].wl,np.array((self.l,Wg[qi,:],Wc[qi,:])).T)

            if store:
                # store
                self.bhe_c[q] = bhe_c[q]
                self.Ag[q] = Ag
                self.Ac[q] = Ac
                if q=='MV':
                    for qi, qq in enumerate(self.qMV): 
                        self.Wg[qq] = Wg[qi,:]
                        self.Wc[qq] = Wc[qi,:]
                        
            if val_return:
                Ags[q] = Ag
                Acs[q] = Ac

        if val_return:
            return Ags, Acs
                

    def loadnorm(self):

        self.Ag, self.Ac, self.Wg, self.Wc = {}, {}, {}, {}

        # load normalization
        for q in self.qlist:
            self.Ag[q], self.Ac[q] = (np.loadtxt(self.f[q].al,usecols=(1,2))).T

        # load optimal weights
        if 'MV' in self.qlist and self.qtype=='lens':
            for q in ['TT','TE','EE','TB','EB']:
                self.Wg[q], self.Wc[q] = (np.loadtxt(self.f[q].wl,usecols=(1,2))).T

        # check BHE responses
        if self.bhe_do and ( not hasattr(self,'bhe_c') or not self.bhe_c ): # true if BHE corrections are not pre-computed
            self.al(output=False)
        


    def qrec(self,qout=None,gtype='k'):
        '''
        Return quadratic estimators
        '''

        Lmax  = self.olmax
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        lcl   = self.lcl

        if qout == None:  qout = self

        # load normalization and weights including BHE coefficients
        self.loadnorm()
        
        # loop for realizations
        if len(self.rlz) <= 0: 
            print('nothing to do for qrec')

        for i in tqdm.tqdm(self.rlz,ncols=100,desc='reconstruction:'):
            
            gmv, cmv = 0., 0.

            # check file exits
            qlist = []
            for q in self.qlist:
                if misctools.check_path(qout.f[q].alm[i],overwrite=self.overwrite,verbose=self.verbose): continue
                qlist.append(q)
            if qlist==[]: continue

            # load cmb alms
            alm = { cmb: self.Fl[cmb][:,None] * pickle.load(open(self.falm[cmb][i],"rb"))[:rlmax+1,:rlmax+1] for cmb in self.mtype }

            for q in tqdm.tqdm(qlist,ncols=100,desc='each quad-comb:',leave=False):

                llm, tlm, slm = 0., 0., 0.
                
                if self.qtype=='lens' or 'lens' in self.bhe_list:
                    if q=='TT':  glm, clm = curvedsky.rec_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],alm['T'],alm['T'],gtype=gtype,nside_t=nside)
                    if q=='TE':  glm, clm = curvedsky.rec_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],alm['T'],alm['E'],gtype=gtype,nside_t=nside)
                    if q=='TB':  glm, clm = curvedsky.rec_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm['T'],alm['B'],gtype=gtype,nside_t=nside)
                    if q=='EE':  glm, clm = curvedsky.rec_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],alm['E'],alm['E'],gtype=gtype,nside_t=nside)
                    if q=='EB':  glm, clm = curvedsky.rec_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm['E'],alm['B'],gtype=gtype,nside_t=nside)
                    if q=='MV':  glm, clm = gmv.copy(), cmv.copy()
                    llm = glm.copy()

                if self.qtype=='rot':
                    if q=='TB':  glm = curvedsky.rec_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm['T'],alm['B'],nside_t=nside)
                    if q=='EB':  glm = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm['E'],alm['B'],nside_t=nside)
                    clm = glm*0.

                if self.qtype=='tau' or 'tau' in self.bhe_list:
                    if q=='TT':  glm = curvedsky.rec_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],alm['T'],alm['T'],nside_t=nside)
                    if q=='EB':  glm = curvedsky.rec_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],alm['E'],alm['B'],nside_t=nside)
                    clm = glm*0.
                    tlm = glm.copy()

                if self.qtype=='src' or 'src' in self.bhe_list:
                    if q=='TT':  glm = curvedsky.rec_src.qtt(Lmax,rlmin,rlmax,alm['T'],alm['T'],nside_t=nside)
                    clm = glm*0.
                    slm = glm.copy()

                if self.qtype=='ilens' or 'ilens' in self.bhe_list:
                    if q=='TE':  glm, clm = curvedsky.rec_ilens.qte(Lmax,rlmin,rlmax,lcl[3,:],alm['T'],alm['E'],gtype=gtype,nside_t=nside)
                    if q=='TB':  glm, clm = curvedsky.rec_ilens.qtb(Lmax,rlmin,rlmax,lcl[3,:],alm['T'],alm['B'],gtype=gtype,nside_t=nside)
                    if q=='EE':  glm, clm = curvedsky.rec_ilens.qee(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],alm['E'],alm['E'],gtype=gtype,nside_t=nside)
                    if q=='EB':  glm, clm = curvedsky.rec_ilens.qeb(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],alm['E'],alm['B'],gtype=gtype,nside_t=nside)
                    if q=='BB':  glm, clm = curvedsky.rec_ilens.qbb(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],alm['B'],alm['B'],gtype=gtype,nside_t=nside)
                    if q=='MV':  glm, clm = gmv.copy(), cmv.copy()
                    llm = glm.copy()

                # normalization correction
                glm *= self.Ag[q][:,None]
                clm *= self.Ac[q][:,None]                

                # Bias hardened estimator
                if self.bhe_do:
                    glm = self.bhe_c[q]['tau'][:,None]*tlm + self.bhe_c[q]['lens'][:,None]*llm + self.bhe_c[q]['src'][:,None]*slm
                
                # save
                pickle.dump((glm,clm),open(qout.f[q].alm[i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

                # MV
                if q in self.qMV and 'MV' in self.qlist:
                    gmv += self.Wg[q][:,None]*glm
                    cmv += self.Wc[q][:,None]*clm

                    
    def n0(self,falm1=None):
        '''
        The N0 bias calculation
        falm1: use different simulation for id1
        '''

        for q in self.qlist:
            if misctools.check_path(self.f[q].n0bs,overwrite=self.overwrite,verbose=self.verbose): return

        # load normalization and weights including BHE coefficients
        self.loadnorm()

        falm0  = self.falm
        Lmax   = self.olmax

        if falm1 is None: falm1 = self.falm

        # power spectrum
        cl = {q: np.zeros((2,Lmax+1)) for q in self.qlist}

        # loop for realizations
        for i in tqdm.tqdm(self.n0rlz,ncols=100,desc='N0 bias:'):

            id0, id1 = 2*i-1, 2*i
            gmv, cmv = 0., 0.

            alm0 = { m: self.Fl[m][:,None]*pickle.load(open(falm0[m][id0],"rb"))[:self.rlmax+1,:self.rlmax+1] for m in self.mtype }
            alm1 = { m: self.Fl[m][:,None]*pickle.load(open(falm1[m][id1],"rb"))[:self.rlmax+1,:self.rlmax+1] for m in self.mtype }

            for q in self.qlist:
            #for q in tqdm.tqdm(self.qlist,ncols=100,desc='each quad-comb:',leave=False):

                if q == 'MV':
                    glm, clm = gmv.copy(), cmv.copy()
                else:
                    X, Y = q[0], q[1]
                    glm, clm = self.qXY(q,Lmax,alm0[X],alm1[X],alm0[Y],alm1[Y])

                if not self.bhe_do:
                    glm *= self.Ag[q][:,None]
                clm *= self.Ac[q][:,None]

                cl[q][0,:] += curvedsky.utils.alm2cl(glm)/(2*self.wn[4]*self.n0sim)
                cl[q][1,:] += curvedsky.utils.alm2cl(clm)/(2*self.wn[4]*self.n0sim)

                # sum for MV, later used in q==MV
                if q in self.qMV and 'MV' in self.qlist:
                    gmv += self.Wg[q][:,None]*glm
                    cmv += self.Wc[q][:,None]*clm

            if self.output_n0bs_rlz: # Output each realization pair for the N0 bias calculation
                cl_rlz = {q: np.zeros((2,Lmax+1)) for q in self.qlist}
                cl_rlz[q][0,:] = curvedsky.utils.alm2cl(glm)/(2*self.wn[4]*self.n0sim)
                cl_rlz[q][1,:] = curvedsky.utils.alm2cl(clm)/(2*self.wn[4]*self.n0sim)
                for q in self.qlist:
                    np.savetxt(self.f[q].n0bs_rlz[i],np.concatenate((self.l[None,:],cl[q])).T)

        for q in self.qlist:
            np.savetxt(self.f[q].n0bs,np.concatenate((self.l[None,:],cl[q])).T)


    def diag_rdn0(self,frcl=None,scl=None,dsmix_only=False,rlz=None):
        '''
        Compute the diagonal-RDN0 for data realization:
          N^diag[Cl^d,Cl^s]
        This routine output the above quantity for each Cl^d. 
        See Hanson et al. (2011) arXiv:1008.4403 for the definition of the diagonal-RDN0. 
        
          scl:        observed cl obtained from simulation average, Cl^s
          frcl:       filenames for data cl
          dsmix_only: extracting data-sim mixed term only
          rlz:        compute only "rlz"-th realization
        
        '''

        self.loadnorm()

        ocl = self.ocl
        if scl is None:
            scl = self.ocl  # use analytic observed cl

        if rlz is None:
            rlzs = self.rlz
        else:
            rlzs = rlz
            
        for i in tqdm.tqdm(rlzs,ncols=100,desc='Diag-RDN0:'):

            if frcl is not None:
                rcl = np.loadtxt(frcl[i],unpack=True,usecols=(1,2,3,4))  # cmb aps for real data
                
            else:
                # load data cmb alms
                rcl = self.__compute_cmb_aps__(i)
        
            rcl[np.where(rcl==0)] = 1e30 # a large number

            # data x data
            cl = ocl**2/rcl
            Ags0, Acs0 = self.al(ocls=cl,output=False,store=False,val_return=True)

            # (data-sim) x (data-sim)
            cl = ocl**2/(scl-rcl)
            Ags1, Acs1 = self.al(ocls=cl,output=False,store=False,val_return=True)

            if dsmix_only:
                scl[np.where(scl==0)] = 1e30 # a large number
                # sim x sim
                cl = ocl**2/scl
                Ags2, Acs2 = self.al(ocls=cl,output=False,store=False,val_return=True)

            for q in self.qlist:

                Ags0[q][np.where(Ags0[q]==0)] = 1e30
                Ags1[q][np.where(Ags1[q]==0)] = 1e30
                Acs0[q][np.where(Acs0[q]==0)] = 1e30
                Acs1[q][np.where(Acs1[q]==0)] = 1e30

                if dsmix_only:
                    Ags2[q][np.where(Ags2[q]==0)] = 1e30
                    Acs2[q][np.where(Acs2[q]==0)] = 1e30
                    n0g = self.Ag[q]**2*(1./Ags0[q]-1./Ags1[q]+1./Ags2[q])
                    n0c = self.Ac[q]**2*(1./Acs0[q]-1./Acs1[q]+1./Acs2[q])
                else:
                    n0g = self.Ag[q]**2*(1./Ags0[q]-1./Ags1[q])
                    n0c = self.Ac[q]**2*(1./Acs0[q]-1./Acs1[q])
                
                np.savetxt(self.f[q].drdn0[i],np.array((self.l,n0g,n0c)).T)

    
    def diag_rdn0_rlz(self,i,ocl=None):
        '''
        Compute diagonal-RDN0 for a specific data realization, i, over a loop for sim cls: 
          N^diag[Cl^d,Cl^s] = <N^diag[Cl^d,Cl^I]>_I
        This routine output N^diag[Cl^d,Cl^I]
        '''

        self.loadnorm()

        if ocl is None:
            ocl = self.ocl

        # compute data aps
        rcl = self.__compute_cmb_aps__(i)
        rcl[np.where(rcl==0)] = 1e30 # a large number

        for I in tqdm.tqdm(self.rlz,ncols=100,desc='Diag-RDN0_rlz:'):

            if I==i: continue

            # prepare sim CMB cls at I-th realization
            scl = self.__compute_cmb_aps__(I)
        
            # data x data
            cl = ocl**2/rcl
            Ags0, Acs0 = self.al(ocls=cl,output=False,store=False,val_return=True)

            # (data-sim) x (data-sim)
            cl = ocl**2/(scl-rcl)
            Ags1, Acs1 = self.al(ocls=cl,output=False,store=False,val_return=True)

            for q in self.qlist:

                Ags0[q][np.where(Ags0[q]==0)] = 1e30
                Ags1[q][np.where(Ags1[q]==0)] = 1e30
                Acs0[q][np.where(Acs0[q]==0)] = 1e30
                Acs1[q][np.where(Acs1[q]==0)] = 1e30

                n0g = self.Ag[q]**2*(1./Ags0[q]-1./Ags1[q])
                n0c = self.Ac[q]**2*(1./Acs0[q]-1./Acs1[q])
                
                np.savetxt(self.f[q].drdn0[i].replace('.dat','_sim'+str(I)+'.dat'),np.array((self.l,n0g,n0c)).T)


    def __compute_cmb_aps__(self,i):
        '''
        Utility function of this object to prepare a set of CMB cls
        '''

        # load cmb alms
        alm = { cmb: pickle.load(open(self.falm[cmb][i],"rb"))[:self.rlmax+1,:self.rlmax+1] for cmb in self.mtype }
        
        # compute aps
        scl = np.zeros((4,self.rlmax+1))
        if 'TT' in self.qlist or 'TE' in self.qlist or 'TB' in self.qlist: 
            scl[0,:] = curvedsky.utils.alm2cl(alm['T'])
        if 'TE' in self.qlist or 'EE' in self.qlist or 'EB' in self.qlist: 
            scl[1,:] = curvedsky.utils.alm2cl(alm['E'])
        if 'TB' in self.qlist or 'EB' in self.qlist or 'BB' in self.qlist: 
            scl[2,:] = curvedsky.utils.alm2cl(alm['B'])
        if 'TE' in self.qlist: 
            scl[3,:] = curvedsky.utils.alm2cl(alm['T'],alm['E'])

        return scl


    
    def rdn0(self,qout=None,falms=None):
        '''
        The full RDN0 bias calculation
        See Namikawa et al. (2013) for the definition of the full RDN0. 
        '''

        # load normalization and weights
        self.loadnorm()

        falm  = self.falm
        Lmax  = self.olmax
        qlist = self.qlist

        if falms is None: falms = self.falm
        if qout  is None:  qout = self

        # load N0
        N0 = { q: np.loadtxt(self.f[q].n0bs,unpack=True,usecols=(1,2)) for q in qlist }

        # compute RDN0
        for i in tqdm.tqdm(self.rlz,ncols=100,desc='RDN0:'):

            # skip RDN0 for sim
            if not self.rd4sim and i!=0: 
                tqdm.tqdm.write('skip rdn0 for sim')
                break

            # avoid overwriting
            Qlist = []
            for q in qlist:
                if misctools.check_path(qout.f[q].rdn0[i],overwrite=self.overwrite,verbose=False):
                    if self.verbose: tqdm.tqdm.write(qout.f[q].rdn0[i]+' exist and is not overwritten')
                    Qlist.append(q)
            if Qlist != []: 
                continue

            # power spectrum
            cl = { q: np.zeros((2,Lmax+1)) for q in qlist }

            # load alm
            almr = { m: self.Fl[m][:,None]*pickle.load(open(falm[m][i],"rb"))[:self.rlmax+1,:self.rlmax+1] for m in self.mtype }

            # loop for sim index (I)
            for I in tqdm.tqdm(range(self.rdmin,self.rdmax+1),ncols=100,desc='inside loop:',leave=False):

                if I==i: continue

                gmv, cmv = 0., 0.

                # load alm
                alms = { m: self.Fl[m][:,None]*pickle.load(open(falms[m][I],"rb"))[:self.rlmax+1,:self.rlmax+1] for m in self.mtype }

                for q in qlist:

                    X, Y = q[0], q[1]

                    if q=='MV':
                        glm, clm = gmv.copy(), cmv.copy()
                    else:
                        glm, clm = self.qXY(q,Lmax,almr[X],alms[X],almr[Y],alms[Y])

                    if not self.bhe_do:
                        glm *= self.Ag[q][:,None]
                    clm *= self.Ac[q][:,None]

                    # add RDN0 power spectrum from this sim index (I)
                    cl[q][0,:] += curvedsky.utils.alm2cl(glm)
                    cl[q][1,:] += curvedsky.utils.alm2cl(clm)

                    # MV
                    if q in self.qMV and 'MV' in qlist:
                        gmv += self.Wg[q][:,None] * glm
                        cmv += self.Wc[q][:,None] * clm

            if self.rdsim>0:
                
                sn = self.rdsim
                if self.rdmin<=i and i<=self.rdmax:  sn = self.rdsim-1
                
                for q in qlist:
                    cl[q] = cl[q]/(self.wn[4]*sn) - N0[q]
                    np.savetxt(qout.f[q].rdn0[i],np.concatenate((self.l[None,:],cl[q])).T)

    
    def rdn0_rlz(self,i,qout=None,falms=None,dsmix_only=False):
        '''
        The full RDN0 bias calculation for each sim realization. 
        '''

        # load normalization and weights
        self.loadnorm()

        falm  = self.falm
        Lmax  = self.olmax
        qlist = self.qlist

        if falms is None: falms = self.falm
        if qout  is None:  qout = self

        # load N0
        if not dsmix_only:
            N0 = { q: np.loadtxt(self.f[q].n0bs,unpack=True,usecols=(1,2)) for q in qlist }

        # load alm
        almr = { m: self.Fl[m][:,None]*pickle.load(open(falm[m][i],"rb"))[:self.rlmax+1,:self.rlmax+1] for m in self.mtype }

        # loop for I
        for I in tqdm.tqdm(range(self.rdmin,self.rdmax+1),ncols=100,desc='RDN0 inside loop:'):

            if I==i: continue
            gmv, cmv = 0., 0.

            # avoid overwriting
            Qlist = []
            for q in qlist:
                filename = qout.f[q].rdn0_rlz[i][I]
                if misctools.check_path(filename,overwrite=self.overwrite,verbose=False):
                    if self.verbose: tqdm.tqdm.write(filename+' exist and is not overwritten')
                    Qlist.append(q)
            if Qlist != []: 
                continue

            # load alm
            alms = { m: self.Fl[m][:,None]*pickle.load(open(falms[m][I],"rb"))[:self.rlmax+1,:self.rlmax+1] for m in self.mtype }

            # RDN0 power spectrum
            cl = { q: np.zeros((2,Lmax+1)) for q in qlist }
    
            for q in qlist:

                X, Y = q[0], q[1]

                if q=='MV':
                    glm, clm = gmv.copy(), cmv.copy()
                else:
                    glm, clm = self.qXY(q,Lmax,almr[X],alms[X],almr[Y],alms[Y])

                if not self.bhe_do:
                    glm *= self.Ag[q][:,None]
                clm *= self.Ac[q][:,None]

                cl[q][0,:] = curvedsky.utils.alm2cl(glm)
                cl[q][1,:] = curvedsky.utils.alm2cl(clm)

                # MV
                if q in self.qMV and 'MV' in qlist:
                    gmv += self.Wg[q][:,None] * glm
                    cmv += self.Wc[q][:,None] * clm

            for q in qlist:
                if not dsmix_only:
                    cl[q] = cl[q]/(self.wn[4]) - N0[q]
                else:
                    cl[q] = cl[q]/(self.wn[4])
                np.savetxt(qout.f[q].rdn0_rlz[i][I],np.concatenate((self.l[None,:],cl[q])).T)

    
    def qXY(self,q,Lmax,Xlm1,Xlm2,Ylm1,Ylm2,gtype='k'):# for N0 and RDN0 estimates
    
        rlmin = self.rlmin
        rlmax = self.rlmax
        nside = self.nside
        lcl   = self.lcl

        glm1, glm2, tlm1, tlm2, slm1, slm2 = 0., 0., 0., 0., 0., 0.
        if self.qtype=='lens' or 'lens' in self.bhe_list:
            if q=='TT':
                glm1, clm1 = curvedsky.rec_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_lens.qtt(Lmax,rlmin,rlmax,lcl[0,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='TE':
                glm1, clm1 = curvedsky.rec_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_lens.qte(Lmax,rlmin,rlmax,lcl[3,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='TB':
                glm1, clm1 = curvedsky.rec_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_lens.qtb(Lmax,rlmin,rlmax,lcl[3,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='EE':
                glm1, clm1 = curvedsky.rec_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_lens.qee(Lmax,rlmin,rlmax,lcl[1,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='EB':
                glm1, clm1 = curvedsky.rec_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_lens.qeb(Lmax,rlmin,rlmax,lcl[1,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)

            if not self.bhe_do:
                return glm1+glm2, clm1+clm2

        if self.qtype=='rot':
            if q=='TB':
                rlm1 = curvedsky.rec_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],Xlm1,Ylm2,nside_t=nside)
                rlm2 = curvedsky.rec_rot.qtb(Lmax,rlmin,rlmax,lcl[3,:],Xlm2,Ylm1,nside_t=nside)
            if q=='EB':
                rlm1 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],Xlm1,Ylm2,nside_t=nside)
                rlm2 = curvedsky.rec_rot.qeb(Lmax,rlmin,rlmax,lcl[1,:],Xlm2,Ylm1,nside_t=nside)

            return rlm1+rlm2, (rlm1+rlm2)*0.

        if self.qtype=='tau' or 'tau' in self.bhe_list:
            if q=='TT':
                tlm1 = curvedsky.rec_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],Xlm1,Ylm2,nside_t=nside)
                tlm2 = curvedsky.rec_tau.qtt(Lmax,rlmin,rlmax,lcl[0,:],Xlm2,Ylm1,nside_t=nside)
            if q=='EB':
                tlm1 = curvedsky.rec_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],Xlm1,Ylm2,nside_t=nside)
                tlm2 = curvedsky.rec_tau.qeb(Lmax,rlmin,rlmax,lcl[1,:],Xlm2,Ylm1,nside_t=nside)

            if not self.bhe_do:
                return tlm1+tlm2, (tlm1+tlm2)*0.

        if self.qtype=='src' or 'src' in self.bhe_list:
            if q=='TT':
                slm1 = curvedsky.rec_src.qtt(Lmax,rlmin,rlmax,Xlm1,Ylm2,nside_t=nside)
                slm2 = curvedsky.rec_src.qtt(Lmax,rlmin,rlmax,Xlm2,Ylm1,nside_t=nside)

            if not self.bhe_do:
                return slm1+slm2, (slm1+slm2)*0.

        if self.qtype=='ilens' or 'ilens' in self.bhe_list:
            if q=='TE':
                glm1, clm1 = curvedsky.rec_ilens.qte(Lmax,rlmin,rlmax,lcl[3,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_ilens.qte(Lmax,rlmin,rlmax,lcl[3,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='TB':
                glm1, clm1 = curvedsky.rec_ilens.qtb(Lmax,rlmin,rlmax,lcl[3,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_ilens.qtb(Lmax,rlmin,rlmax,lcl[3,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='EE':
                glm1, clm1 = curvedsky.rec_ilens.qee(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_ilens.qee(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='EB':
                glm1, clm1 = curvedsky.rec_ilens.qeb(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_ilens.qeb(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)
            if q=='BB':
                glm1, clm1 = curvedsky.rec_ilens.qbb(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],Xlm1,Ylm2,gtype=gtype,nside_t=nside)
                glm2, clm2 = curvedsky.rec_ilens.qbb(Lmax,rlmin,rlmax,lcl[1,:]-lcl[2,:],Xlm2,Ylm1,gtype=gtype,nside_t=nside)

            if not self.bhe_do:
                return glm1+glm2, clm1+clm2


        # Bias hardened estimator (This alm is already normalized)
        if self.bhe_do:
            alm1 = self.bhe_c[q]['tau'][:,None]*tlm1 + self.bhe_c[q]['lens'][:,None]*glm1 + self.bhe_c[q]['src'][:,None]*slm1
            alm2 = self.bhe_c[q]['tau'][:,None]*tlm2 + self.bhe_c[q]['lens'][:,None]*glm2 + self.bhe_c[q]['src'][:,None]*slm2
            return alm1+alm2, (alm1+alm2)*0.


    def mean(self):

        Lmax = self.olmax
        rlz  = np.linspace(self.mfmin,self.mfmax,self.mfmax-self.mfmin+1,dtype=np.int32)

        for q in self.qlist:

            if misctools.check_path(self.f[q].MFalm,overwrite=self.overwrite,verbose=self.verbose): continue

            mfg, mfc = 0., 0.
            for I in tqdm.tqdm(rlz,ncols=100,desc='mean-field (all sim average): ('+q+')'):
                mfgi, mfci = pickle.load(open(self.f[q].alm[I],"rb"))
                mfg += mfgi/self.mfsim
                mfc += mfci/self.mfsim

            pickle.dump((mfg,mfc),open(self.f[q].MFalm,"wb"),protocol=pickle.HIGHEST_PROTOCOL)

            # compute mf cls
            if self.verbose:  print('cl for mean field bias')
            cl = np.zeros((2,Lmax+1))
            cl[0,:] = curvedsky.utils.alm2cl(mfg)/self.wn[4]
            cl[1,:] = curvedsky.utils.alm2cl(mfc)/self.wn[4]
            np.savetxt(self.f[q].MFcl,np.concatenate((self.l[None,:],cl)).T)


    def mean_rlz(self):

        Lmax = self.olmax

        for q in self.qlist:

            # counting missing files
            filen = 0
            for i in self.rlz:
                if misctools.check_path(self.f[q].mfalm[i],overwrite=self.overwrite,verbose=self.verbose): continue
                filen += 1
            if filen == 0: 
                continue  # do nothing below

            mglm = np.zeros((Lmax+1,Lmax+1),dtype=np.complex)
            mclm = np.zeros((Lmax+1,Lmax+1),dtype=np.complex)

            for I in tqdm.tqdm(self.mfrlz,ncols=100,desc='mean-field: load reconstructed alms ('+q+')'):
                glm, clm = pickle.load(open(self.f[q].alm[I],"rb"))
                mglm += glm/self.mfsim
                mclm += clm/self.mfsim

            for i in tqdm.tqdm(self.rlz,ncols=100,desc='mean-field: alm for each rlz ('+q+')'):

                if misctools.check_path(self.f[q].mfalm[i],overwrite=self.overwrite,verbose=self.verbose): continue
        
                if i>=self.mfmin and i<=self.mfmax:
                    glm, clm = pickle.load(open(self.f[q].alm[i],"rb"))
                    mfg = mglm - glm/self.mfsim
                    mfc = mclm - clm/self.mfsim
                    mfg *= self.mfsim/(self.mfsim-1.)
                    mfc *= self.mfsim/(self.mfsim-1.)
                else:
                    mfg = 1.*mglm
                    mfc = 1.*mclm

                pickle.dump((mfg,mfc),open(self.f[q].mfalm[i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

            # compute mf cls
            for i in tqdm.tqdm(self.rlz,ncols=100,desc='mean-field: aps for each rlz ('+q+')'):

                if misctools.check_path(self.f[q].mfcl[i],overwrite=self.overwrite,verbose=self.verbose): continue

                mfg, mfc = pickle.load(open(self.f[q].mfalm[i],"rb"))

                cl = np.zeros((2,Lmax+1))
                cl[0,:] = curvedsky.utils.alm2cl(mfg) / self.wn[4]
                cl[1,:] = curvedsky.utils.alm2cl(mfc) / self.wn[4]
                np.savetxt(self.f[q].mfcl[i],np.concatenate((self.l[None,:],cl)).T)



    def qrec_flow(self,run=[],qout=None,falm_n0bs=None,falm_sim=None):
        '''
        falm_n0:  this different simulation (falm1) will be used for id1 in N0
        falm_sim: this alms will be used simulation alms in RDN0
        '''

        # create directory to store intermediate files
        if 'directory' in run:
            self.create_directory()

        # set filtering
        if run:
            self.cinvfilter()

        # normalization
        if 'norm' in run:
            self.al()

        # quadratic estimators
        if 'qrec' in run:
            self.qrec()

        # Realization-independent N0
        if 'n0' in run:
            self.n0(falm1=falm_n0bs)

        # Realization-dependent N0
        if 'rdn0' in run:
            self.rdn0(qout=qout,falms=falm_sim)

        # mean-field bias
        if 'mean' in run:
            self.mean()

        # mean-field bias
        if 'mean_rlz' in run:
            self.mean_rlz()




# Define quadratic estimator names
class quad_fname_cross():

    def __init__(self,qobj,qest):

        # qtype is the type of mode coupling, such as lens, rot, etc
        qrdn = qobj.qrdn
        qaps = qobj.qaps

        qtag = qest + '_' + qobj.qobj0.cmbtag + qobj.qobj1.cmbtag + qobj.qobj0.ltag + qobj.qobj1.ltag + qobj.qobj0.qtagext + qobj.qobj1.qtagext

        # N0 bias
        self.n0bs = qaps+'n0_'+qtag+'_n'+str(qobj.n0sim).zfill(3)+'.dat'

        # reconstructed alm and RDN0
        self.rdn0 = [qrdn+'rdn0_'+qtag+'_n'+str(qobj.rdsim).zfill(3)+'_'+x+'.dat' for x in qobj.ids]

        # RDN0 for each sim rlz
        self.rdn0_rlz = [ {I: qrdn+'rdn0_'+qtag+'_rd'+str(I).zfill(3)+'_'+x+'.dat' for I in qobj.rdids } for x in qobj.ids ]

        # diagonal RDN0
        self.drdn0 = [qrdn+'drdn0_'+qtag+'_'+x+'.dat' for x in qobj.ids]


        
class quad_cross(): # for phi cross-spectrum between two different CMB data

    def __init__(self, qobj0, qobj1, 
                 ids=[],
                 olmax=2048, nside=2048,
                 n0min=1, n0max=50, rdmin=1, rdmax=100, 
                 qlist=None, qtype='', wn=None, bhe=None, rd4sim=False, 
                 overwrite=False, verbose=True, 
                 root=''
                ):
        
        #//// get parameters ////#
        conf = misctools.load_config('QUADREC')
        self.ids = ids

        self.qobj0 = qobj0
        self.qobj1 = qobj1
        
        self.olmax = conf.getint('olmax',olmax)

        # start, stop rlz of N0 bias
        self.n0min  = conf.getint('n0min',n0min)
        self.n0max  = conf.getint('n0max',n0max)
        self.n0sim  = self.n0max - self.n0min + 1
        self.n0rlz  = np.linspace(self.n0min,self.n0max,self.n0sim,dtype=np.int32)

        # start, stop rlz of RDN0 bias
        self.rdmin  = conf.getint('rdmin',rdmin)
        self.rdmax  = conf.getint('rdmax',rdmax)
        self.rdsim  = self.rdmax - self.rdmin + 1
        self.rdids  = np.arange(self.rdmin,self.rdmax+1)  # index for RDN0 sim
        self.rd4sim = conf.getboolean('rd4sim',rd4sim)  # whether RD calculation for sim

        self.qlist = qlist
        self.qtype = qtype
        self.mtype = set_mtype(self.qlist)
        
        # window normalization correction
        if wn is None:
            self.wn = np.ones(5)
        else:
            self.wn = wn

        #//// Bias Herdened Estimators ////
        setup_bhe(self,bhe)

        #//// File and Directory ////#
        self.root = root
        self.qrdn = self.root + self.qtype + '/rdn0/'    # RDN0
        self.qaps = self.root + self.qtype + '/aps/'     # Cls of reconstructed alm

        #setup filename
        self.f = {q: quad_fname_cross(self,q) for q in ['TT','TE','TB','EE','EB','BB','MV']}

        #//// Misc ////#
        self.overwrite = overwrite
        self.verbose = verbose
 

        
    def n0(self):
        '''
        Cross qobj0 x qobj1 
        N0 for X^{A,1}Y^{A,2} Z^{B,1}W^{B,2} + X^{A,1}Y^{A,2} Z^{B,2}W^{B,1} 
        '''
        
        qobj0 = self.qobj0
        qobj1 = self.qobj1

        # load normalization and weights
        qobj0.loadnorm()
        qobj1.loadnorm()

        # power spectrum
        cl = {q: np.zeros((2,self.olmax+1)) for q in self.qlist}

        # loop for realizations
        for i in tqdm.tqdm(self.n0rlz,ncols=100,desc='N0 bias (cross):'):

            id0, id1 = 2*i-1, 2*i
            gmv, cmv = 0., 0.

            alm0 = { m: qobj0.Fl[m][:,None] * pickle.load(open(qobj0.falm[m][id0],"rb"))[:qobj0.rlmax+1,:qobj0.rlmax+1] for m in self.mtype }
            alm1 = { m: qobj0.Fl[m][:,None] * pickle.load(open(qobj0.falm[m][id1],"rb"))[:qobj0.rlmax+1,:qobj0.rlmax+1] for m in self.mtype }
            blm0 = { m: qobj1.Fl[m][:,None] * pickle.load(open(qobj1.falm[m][id0],"rb"))[:qobj1.rlmax+1,:qobj1.rlmax+1] for m in self.mtype }
            blm1 = { m: qobj1.Fl[m][:,None] * pickle.load(open(qobj1.falm[m][id1],"rb"))[:qobj1.rlmax+1,:qobj1.rlmax+1] for m in self.mtype }

            for q in tqdm.tqdm(self.qlist,ncols=100,desc='each quad-comb:',leave=False):

                glm0, glm1, clm0, clm1 = self.qXY(q,Lmax,alm0,alm1,blm0,blm1)

                if not qobj0.bhe_do: glm0 *= qobj0.Ag[q][:,None]
                if not qobj1.bhe_do: glm1 *= qobj1.Ag[q][:,None]
                clm0 *= qobj0.Ac[q][:,None]
                clm1 *= qobj1.Ac[q][:,None]

                cl[q][0,:] += curvedsky.utils.alm2cl(glm0,glm1)/(2*self.wn[4]*self.n0sim)
                cl[q][1,:] += curvedsky.utils.alm2cl(clm0,clm1)/(2*self.wn[4]*self.n0sim)


        for q in self.qlist:
            print ('save N0 data')
            oL = np.linspace(0,self.olmax,self.olmax+1)
            np.savetxt(self.f[q].n0bs,np.concatenate((oL[None,:],cl[q])).T)



    def rdn0(self):
        '''
        The sim-data-mixed term of the RDN0 bias calculation
        '''

        qobj0 = self.qobj0
        qobj1 = self.qobj1

        # load normalization and weights
        qobj0.loadnorm()
        qobj1.loadnorm()

        # maximum multipole of output
        Lmax  = self.olmax
        oL = np.linspace(0,Lmax,Lmax+1)

        # load N0
        N0 = {q: np.loadtxt(self.f[q].n0bs,unpack=True,usecols=(1,2)) for q in self.qlist}

        # compute RDN0
        for i in tqdm.tqdm(self.rdids,ncols=100,desc='RDN0 (Cross):'):

            # power spectrum
            cl = {q: np.zeros((2,Lmax+1)) for q in self.qlist}

            # load alm
            almr = { m: qobj0.Fl[m][:,None]*pickle.load(open(qobj0.falm[m][i],"rb"))[:qobj0.rlmax+1,:qobj0.rlmax+1] for m in self.mtype}
            blmr = { m: qobj1.Fl[m][:,None]*pickle.load(open(qobj1.falm[m][i],"rb"))[:qobj1.rlmax+1,:qobj1.rlmax+1] for m in self.mtype}

            # loop for sim
            for I in tqdm.tqdm(range(self.rdmin,self.rdmax+1),ncols=100,desc='inside loop:',leave=False):

                # load alm
                alms = { m: qobj0.Fl[m][:,None]*pickle.load(open(qobj0.falm[m][I],"rb"))[:qobj0.rlmax+1,:qobj0.rlmax+1] for m in self.mtype }
                blms = { m: qobj1.Fl[m][:,None]*pickle.load(open(qobj1.falm[m][I],"rb"))[:qobj1.rlmax+1,:qobj1.rlmax+1] for m in self.mtype }

                for q in self.qlist:

                    if I==i: continue

                    glm0, glm1, clm0, clm1 = self.qXY(q,Lmax,almr,alms,blmr,blms)

                    if not qobj0.bhe_do: glm0 *= qobj0.Ag[q][:,None]
                    if not qobj1.bhe_do: glm1 *= qobj1.Ag[q][:,None]
                    clm0 *= qobj0.Ac[q][:,None]
                    clm1 *= qobj1.Ac[q][:,None]

                    cl[q][0,:] += curvedsky.utils.alm2cl(glm0,glm1)
                    cl[q][1,:] += curvedsky.utils.alm2cl(clm0,clm1)


            if self.rdsim>0:
                
                sn = self.rdsim
                if self.rdmin<=i and i<=self.rdmax:  sn = self.rdsim-1

                for q in self.qlist:
                    cl[q] = cl[q]/(self.wn[4]*sn) - N0[q]
                    np.savetxt(self.f[q].rdn0[i],np.concatenate((oL[None,:],cl[q])).T)

    
    def rdn0_rlz(self,i,qout=None,falms=None,dsmix_only=False,alpha=0.):

        qobj0 = self.qobj0
        qobj1 = self.qobj1

        # load normalization and weights
        qobj0.loadnorm()
        qobj1.loadnorm()

        # maximum multipole of output
        Lmax  = self.olmax
        oL = np.linspace(0,Lmax,Lmax+1)

        # load N0
        if not dsmix_only:
            N0 = {q: np.loadtxt(self.f[q].n0bs,unpack=True,usecols=(1,2)) for q in self.qlist}

        # power spectrum
        cl = {q: np.zeros((2,Lmax+1)) for q in self.qlist}

        # load alm
        almr = { m: qobj0.Fl[m][:,None]*pickle.load(open(qobj0.falm[m][i],"rb"))[:qobj0.rlmax+1,:qobj0.rlmax+1] for m in self.mtype}
        blmr = { m: qobj1.Fl[m][:,None]*pickle.load(open(qobj1.falm[m][i],"rb"))[:qobj1.rlmax+1,:qobj1.rlmax+1] for m in self.mtype}

        # loop for sim
        for I in tqdm.tqdm(self.rdids,ncols=100,desc='inside loop:',leave=False):

            # load alm
            alms = { m: qobj0.Fl[m][:,None]*pickle.load(open(qobj0.falm[m][I],"rb"))[:qobj0.rlmax+1,:qobj0.rlmax+1] for m in self.mtype }
            blms = { m: qobj1.Fl[m][:,None]*pickle.load(open(qobj1.falm[m][I],"rb"))[:qobj1.rlmax+1,:qobj1.rlmax+1] for m in self.mtype }

            for q in self.qlist:

                if I==i: continue

                glm0, glm1, clm0, clm1 = self.qXY(q,Lmax,almr,alms,blmr,blms,alpha)

                if not qobj0.bhe_do: glm0 *= qobj0.Ag[q][:,None]
                if not qobj1.bhe_do: glm1 *= qobj1.Ag[q][:,None]
                clm0 *= qobj0.Ac[q][:,None]
                clm1 *= qobj1.Ac[q][:,None]

                cl[q][0,:] += curvedsky.utils.alm2cl(glm0,glm1)
                cl[q][1,:] += curvedsky.utils.alm2cl(clm0,clm1)

                if not dsmix_only:
                    cl[q] = cl[q]/(self.wn[4]) - N0[q]
                else:
                    cl[q] = cl[q]/(self.wn[4])
                np.savetxt(self.f[q].rdn0_rlz[i][I],np.concatenate((oL[None,:],cl[q])).T)
                    
    
    def qXY(self,q,Lmax,almr,alms,blmr,blms,alpha=0.,gtype='k'):
        
        # lcl could have different rlmax

        qobj0 = self.qobj0
        qobj1 = self.qobj1
        q1, q2 = q[0], q[1]

        glm01a, glm10a, tlm01a, tlm10a, slm01a, slm10a = 0., 0., 0., 0., 0., 0.
        glm01b, glm10b, tlm01b, tlm10b, slm01b, slm10b = 0., 0., 0., 0., 0., 0.
        
        almr1 = almr[q1]
        almr2 = almr[q2]
        alms1 = alms[q1]
        alms2 = alms[q2]
        blmr1 = blmr[q1]
        blmr2 = blmr[q2]
        blms1 = blms[q1]
        blms2 = blms[q2]

        if alpha != 0:
            alms1 = alms[q1] - alpha*blms[q1]
            alms2 = alms[q2] - alpha*blms[q2]
            blms1 = alms[q1] + alpha*blms[q1]
            blms2 = alms[q2] + alpha*blms[q2]

        if self.qtype=='lens' or 'lens' in self.bhe_list:
            if q=='TT':
                glm01a, clm01a = curvedsky.rec_lens.qtt(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[0,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                glm10a, clm10a = curvedsky.rec_lens.qtt(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[0,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                glm01b, clm01b = curvedsky.rec_lens.qtt(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[0,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                glm10b, clm10b = curvedsky.rec_lens.qtt(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[0,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)
            if q=='TE':
                glm01a, clm01a = curvedsky.rec_lens.qte(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[3,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                glm10a, clm10a = curvedsky.rec_lens.qte(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[3,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                glm01b, clm01b = curvedsky.rec_lens.qte(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[3,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                glm10b, clm10b = curvedsky.rec_lens.qte(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[3,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)
            if q=='TB':
                glm01a, clm01a = curvedsky.rec_lens.qtb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[3,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                glm10a, clm10a = curvedsky.rec_lens.qtb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[3,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                glm01b, clm01b = curvedsky.rec_lens.qtb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[3,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                glm10b, clm10b = curvedsky.rec_lens.qtb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[3,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)
            if q=='EE':
                glm01a, clm01a = curvedsky.rec_lens.qee(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                glm10a, clm10a = curvedsky.rec_lens.qee(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                glm01b, clm01b = curvedsky.rec_lens.qee(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                glm10b, clm10b = curvedsky.rec_lens.qee(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)
            if q=='EB':
                glm01a, clm01a = curvedsky.rec_lens.qeb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                glm10a, clm10a = curvedsky.rec_lens.qeb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                glm01b, clm01b = curvedsky.rec_lens.qeb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                glm10b, clm10b = curvedsky.rec_lens.qeb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)

            if not self.bhe_do:
                return glm01a+glm10a, glm01b+glm10b, clm01a+clm10a, clm01b+clm10b

        if self.qtype=='rot':
            if q=='EB':
                rlm01a = curvedsky.rec_rot.qeb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],almr1,alms2,nside_t=qobj0.nside)
                rlm10a = curvedsky.rec_rot.qeb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],alms1,almr2,nside_t=qobj0.nside)
                rlm01b = curvedsky.rec_rot.qeb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blmr1,blms2,nside_t=qobj1.nside)
                rlm10b = curvedsky.rec_rot.qeb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blms1,blmr2,nside_t=qobj1.nside)

            return rlm01a+rlm10a, rlm01b+rlm10b, (rlm01a+rlm10a)*0., (rlm01b+rlm10b)*0.


        if self.qtype=='tau' or 'tau' in self.bhe_list:
            if q=='TT':
                tlm01a = curvedsky.rec_tau.qtt(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[0,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                tlm10a = curvedsky.rec_tau.qtt(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[0,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                tlm01b = curvedsky.rec_tau.qtt(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[0,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                tlm10b = curvedsky.rec_tau.qtt(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[0,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)
            if q=='EB':
                tlm01a = curvedsky.rec_tau.qeb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                tlm10a = curvedsky.rec_tau.qeb(Lmax,qobj0.rlmin,qobj0.rlmax,qobj0.lcl[1,:],alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                tlm01b = curvedsky.rec_tau.qeb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                tlm10b = curvedsky.rec_tau.qeb(Lmax,qobj1.rlmin,qobj1.rlmax,qobj1.lcl[1,:],blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)

            if not self.bhe_do:
                return tlm01a+tlm10a, tlm01b+tlm10b, (tlm01a+tlm10a)*0., (tlm01b+tlm10b)*0

        if self.qtype=='src' or 'src' in self.bhe_list:
            if q=='TT':
                slm01a = curvedsky.rec_src.qtt(Lmax,qobj0.rlmin,qobj0.rlmax,almr1,alms2,gtype=gtype,nside_t=qobj0.nside)
                slm10a = curvedsky.rec_src.qtt(Lmax,qobj0.rlmin,qobj0.rlmax,alms1,almr2,gtype=gtype,nside_t=qobj0.nside)
                slm01b = curvedsky.rec_src.qtt(Lmax,qobj1.rlmin,qobj1.rlmax,blmr1,blms2,gtype=gtype,nside_t=qobj1.nside)
                slm10b = curvedsky.rec_src.qtt(Lmax,qobj1.rlmin,qobj1.rlmax,blms1,blmr2,gtype=gtype,nside_t=qobj1.nside)

            if not self.bhe_do:
                return slm01a+slm10a, slm01b+slm10b, (slm01a+slm10a)*0., (slm01b+slm10b)*0

        # Bias hardened estimator (This alm is already normalized)
        if self.bhe_do:
            alm01a = qobj0.bhe_c[q]['tau'][:,None]*tlm01a + qobj0.bhe_c[q]['lens'][:,None]*glm01a + qobj0.bhe_c[q]['src'][:,None]*slm01a
            alm10a = qobj0.bhe_c[q]['tau'][:,None]*tlm10a + qobj0.bhe_c[q]['lens'][:,None]*glm10a + qobj0.bhe_c[q]['src'][:,None]*slm10a
            alm01b = qobj1.bhe_c[q]['tau'][:,None]*tlm01b + qobj1.bhe_c[q]['lens'][:,None]*glm01b + qobj1.bhe_c[q]['src'][:,None]*slm01b
            alm10b = qobj1.bhe_c[q]['tau'][:,None]*tlm10b + qobj1.bhe_c[q]['lens'][:,None]*glm10b + qobj1.bhe_c[q]['src'][:,None]*slm10b
            
            return alm01a+alm10a, alm01b+alm10b, (alm01a+alm10a)*0., (alm01b+alm10b)*0.


# ////////////////////// #
# Some useful functions  #
# ////////////////////// #

def reconstruction(droot,ids,rlz=[],getobj=True,run=[],**kwargs):

    """
    Reconstructing CMB lensing potential and its curl mode

    Args:
        :droot (*str*): Root directory of files to be saved
        :ids[] (*str*): 1D array of strings identifying rlz, e.g. 0001, 0002, ...

    Args(optional):
        :rlz[] (*int*): 1D array of integers identifying rlz to be computed
        :stag (*str*): An optional string in order to specify which type of input CMB is used
        :getobj (*bool*): Retrun objects of reconstruction class, containing parameters, filenames
        :run (*str*): Running reconstruction codes
        :**kwargs : Other optional parameters for quad object
    """

    # read parameters
    qobj = quad(rlz=rlz,root=droot,ids=ids,**kwargs)

    # Main calculation
    qobj.qrec_flow(run=run)

    # Return parameters, filenames
    if getobj:
        return qobj

    
def load_rec_alm(qobj,q,rlz,mean_sub=True,mean_rlz=False):
    
    glm, clm = pickle.load(open(qobj.f[q].alm[rlz],"rb"))

    if mean_sub:
        if mean_rlz: # rlz dependent mean field
            mfg, mfc = pickle.load(open(qobj.f[q].mfalm[rlz],"rb"))
            glm -= mfg
            clm -= mfc
        else: # rlz-mean mf
            mfg, mfc = pickle.load(open(qobj.f[q].MFalm,"rb"))
            if rlz==0: 
                glm -= mfg           
                clm -= mfc
            else:
                glm = glm*(1.+1./qobj.mfsim) - mfg           
                clm = clm*(1.+1./qobj.mfsim) - mfc
        
    return glm, clm
            

def cinv_empirical_fltr(lcl,wcl,cnl,ep=1e-30):
    
    # quality factor defined in Planck 2015 lensing paper
    # T' = Q T^f = Q/(cl+nl) * (T+n)/sqrt(Q)
    
    # wcl --- Wiener-filtered cl
    # cnl --- 1D diagonal signal + noise spectrum
    
    if np.ndim(lcl) == 1:
        Ql  = lcl**2/(wcl*cnl+ep**2)
        ocl = np.reshape( cnl/(Ql+ep), (1,len(lcl)) )
        ifl = np.reshape( lcl/(Ql+ep), (1,len(lcl)) )
        
    if np.ndim(lcl) == 2:
        Ql  = (lcl[0:4,:])**2/(wcl*cnl+ep**2)
        ocl = cnl/(Ql+ep) # corrected observed cl  (obs TE is not used)
        ifl = lcl[0:3,:]/(Ql[:3]+ep) # remove theory signal cl in wiener filter
    
    return ocl, ifl

