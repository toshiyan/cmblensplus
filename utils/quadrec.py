
import sys
import numpy as np
import curvedsky
from memory_profiler import profile

if sys.version_info[:3] > (3,0):
  import pickle
elif sys.version_info[:3] > (2,5,2):
  import cPickle as pickle


# compute normalization
@profile
def al(p,f,r):
  '''
  Return normalization of the quadratic estimators
  '''

  for q in p.qlist:

    if p.qtype=='lens':
      if q=='TT': Ag, Ac = curvedsky.norm_lens.qtt(p.lmax,p.rlmin,p.rlmax,r.lcl[0,:],r.oc[0,:],gtype='k')
      if q=='TE': Ag, Ac = curvedsky.norm_lens.qte(p.lmax,p.rlmin,p.rlmax,r.lcl[3,:],r.oc[0,:],r.oc[1,:],gtype='k')
      if q=='EE': Ag, Ac = curvedsky.norm_lens.qee(p.lmax,p.rlmin,p.rlmax,r.lcl[1,:],r.oc[1,:],gtype='k')
      if q=='TB': Ag, Ac = curvedsky.norm_lens.qtb(p.lmax,p.rlmin,p.rlmax,r.lcl[3,:],r.oc[0,:],r.oc[2,:],gtype='k')
      if q=='EB': Ag, Ac = curvedsky.norm_lens.qeb(p.lmax,p.rlmin,p.rlmax,r.lcl[1,:],r.oc[1,:],r.oc[2,:],gtype='k')
      if q=='MV':
        ag, ac, Wg, Wc = curvedsky.norm_lens.qall(p.qDO,p.lmax,p.rlmin,p.rlmax,r.lcl,r.oc,gtype='k')
        Ag, Ac = ag[5,:], ac[5,:]

    if p.qtype=='rot':
      if q=='EB':
        Ac = np.zeros(p.lmax+1)
        Ag = curvedsky.norm_rot.qeb(p.lmax,p.rlmin,p.rlmax,r.lcl[1,:],r.oc[1,:],r.oc[2,:])

    # save
    np.savetxt(f.quad[q].al,np.array((r.eL,Ag,Ac)).T)
    if q=='MV' and p.qtype=='lens': 
      for qi, qq in enumerate(['TT','TE','EE','TB','EB']): np.savetxt(f.quad[qq].wl,np.array((r.eL,Wg[qi,:],Wc[qi,:])).T)



def loadnorm(p,files):
  Ag = {}
  Ac = {}
  Wg = {}
  Wc = {}
  for q in p.qlist:
    Ag[q], Ac[q] = np.loadtxt(files[q].al,unpack=True,usecols=(1,2))
  if 'MV' in p.qlist and p.qtype=='lens': 
    for qi, qq in enumerate(['TT','TE','EE','TB','EB']):  Wg[qq], Wc[qq] = np.loadtxt(files[q].wl,unpack=True,usecols=(1,2))

  return Ag, Ac, Wg, Wc



@profile
def qrec(p,falm,fquad,r):
  '''
  Return quadratic estimators
  '''

  # load normalization and weights
  Ag, Ac, Wg, Wc = loadnorm(p,fquad)

  # loop for realizations
  for i in range(p.snmin,p.snmax):
    print(i)

    gmv = 0.
    cmv = 0.

    for q in p.qlist:

      if 'T' in q:  Talm = r.Fl['T'] * pickle.load(open(falm['T'][i],"rb"))
      if 'E' in q:  Ealm = r.Fl['E'] * pickle.load(open(falm['E'][i],"rb"))
      if 'B' in q:  Balm = r.Fl['B'] * pickle.load(open(falm['B'][i],"rb"))

      if p.qtype=='lens':
        if q=='TT':  glm, clm = curvedsky.rec_lens.qtt(p.lmax,p.rlmin,p.rlmax,r.lcl[0,:],Talm,Talm,gtype='k',nside=p.nsidet)
        if q=='TE':  glm, clm = curvedsky.rec_lens.qte(p.lmax,p.rlmin,p.rlmax,r.lcl[3,:],Talm,Ealm,gtype='k',nside=p.nsidet)
        if q=='TB':  glm, clm = curvedsky.rec_lens.qtb(p.lmax,p.rlmin,p.rlmax,r.lcl[3,:],Talm,Balm,gtype='k',nside=p.nsidet)
        if q=='EE':  glm, clm = curvedsky.rec_lens.qee(p.lmax,p.rlmin,p.rlmax,r.lcl[1,:],Ealm,Ealm,gtype='k',nside=p.nsidet)
        if q=='EB':  glm, clm = curvedsky.rec_lens.qeb(p.lmax,p.rlmin,p.rlmax,r.lcl[1,:],Ealm,Balm,gtype='k',nside=p.nsidet)
        if q=='MV':  glm, clm = gmv, cmv

      if p.qtype=='rot':
        if q=='EB':  glm = curvedsky.rec_rot.qeb(p.lmax,p.rlmin,p.rlmax,r.lcl[1,:],Ealm,Balm,nside=p.nsidet)
        clm = glm*0.

      glm *= Ag[q][:,None]
      clm *= Ac[q][:,None]
      pickle.dump((glm,clm),open(fquad[q].alm[i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

      # T+P
      if q in p.qMV and 'MV' in p.qlist:
        gmv += Wg[q][:,None]*glm
        cmv += Wc[q][:,None]*clm


@profile
def n0(p,falm,fquad,r):
  '''
  The N0 bias calculation
  '''

  # load normalization and weights
  Ag, Ac, Wg, Wc = loadnorm(p,fquad)

  # power spectrum
  cl ={}
  for q in p.qlist:
    cl[q] = np.zeros((2,p.lmax+1))

  # loop for realizations
  for i in range(p.snn0):
    print (2*i+1, 2*i+2)

    gmv = 0.
    cmv = 0.

    for q in p.qlist:

      if q == 'MV':
        glm, clm = gmv, cmv
      else:
        q1, q2 = q[0], q[1]
        print(q1,q2)
        alm1 = r.Fl[q1] * pickle.load(open(falm[q1][2*i+1],"rb"))
        alm2 = r.Fl[q1] * pickle.load(open(falm[q1][2*i+2],"rb"))
        if q1 == q2:
          blm1 = alm1
          blm2 = alm2
        else:
          blm1 = r.Fl[q2] * pickle.load(open(falm[q2][2*i+1],"rb"))
          blm2 = r.Fl[q2] * pickle.load(open(falm[q2][2*i+2],"rb"))
        glm, clm = qXY(q,p,r.lcl,alm1,alm2,blm1,blm2)

      glm *= Ag[q][:,None]
      clm *= Ac[q][:,None]

      cl[q][0,:] += curvedsky.utils.alm2cl(p.lmax,glm)/(2*r.w4*p.snn0)
      cl[q][1,:] += curvedsky.utils.alm2cl(p.lmax,clm)/(2*r.w4*p.snn0)

      # T+P
      if q in p.qMV and 'MV' in p.qlist:
        gmv += Wg[q][:,None]*glm
        cmv += Wc[q][:,None]*clm


  for q in p.qlist:

    if p.snn0>0:
      print ('save N0 data')
      np.savetxt(fquad[q].n0bl,np.concatenate((r.eL[None,:],cl[q])).T)


@profile
def rdn0(p,falm,fquad,r):
  '''
  The 1st set of the RDN0 bias calculation
  '''

  # load normalization and weights
  Ag, Ac, Wg, Wc = loadnorm(p,fquad)

  # load N0
  N0 = {}
  for q in p.qlist:
    N0[q] = np.loadtxt(fquad[q].n0bl,unpack=True,usecols=(1,2))

  # compute RDN0
  for i in range(p.snmin,p.snmax):
    print(i)

    # power spectrum
    cl = {}
    for q in p.qlist:
      cl[q] = np.zeros((2,p.lmax+1))

    # load alm
    almr = {}
    for cmb in ['T','E','B']:
      almr[cmb] = r.Fl[cmb]*pickle.load(open(falm[cmb][i],"rb"))

    # loop for I
    for I in range(1,p.snrd+1):

      gmv = 0.
      cmv = 0.

      # load alm
      alms = {}
      for cmb in ['T','E','B']:
        alms[cmb] = r.Fl[cmb]*pickle.load(open(falm[cmb][I],"rb"))

      for q in p.qlist:

        q1, q2 = q[0], q[1]

        if I==i: continue
        print(I)

        if q=='MV':
          glm, clm = gmv, cmv
        else:
          glm, clm = qXY(q,p,r.lcl,almr[q1],alms[q1],almr[q2],alms[q2])

        glm *= Ag[q][:,None]
        clm *= Ac[q][:,None]

        cl[q][0,:] += curvedsky.utils.alm2cl(p.lmax,glm)
        cl[q][1,:] += curvedsky.utils.alm2cl(p.lmax,clm)

      # T+P
      if q in p.qMV and 'MV' in p.qlist:
        gmv += Wg[q][:,None]*glm
        cmv += Wc[q][:,None]*clm


  if p.snrd>0:
    if i==0:  sn = p.snrd
    if i!=0:  sn = p.snrd-1
    for q in p.qlist:
      cl[q] = cl[q]/(r.w4*sn) - N0[q]
      print ('save RDN0')
      np.savetxt(fquad[q].rdn0[i],np.concatenate((r.eL[None,:],cl[q])).T)


def qXY(q,p,lcl,alm1,alm2,blm1,blm2):

  if p.qtype=='lens':

    if q=='TT':
      glm1, clm1 = curvedsky.rec_lens.qtt(p.lmax,p.rlmin,p.rlmax,lcl[0,:],alm1,alm2,gtype='k',nside=p.nsidet)
      glm2, clm2 = curvedsky.rec_lens.qtt(p.lmax,p.rlmin,p.rlmax,lcl[0,:],alm2,alm1,gtype='k',nside=p.nsidet)

    if q=='TE':
      glm1, clm1 = curvedsky.rec_lens.qte(p.lmax,p.rlmin,p.rlmax,lcl[3,:],alm1,blm2,gtype='k',nside=p.nsidet)
      glm2, clm2 = curvedsky.rec_lens.qte(p.lmax,p.rlmin,p.rlmax,lcl[3,:],alm2,blm1,gtype='k',nside=p.nsidet)

    if q=='TB':
      glm1, clm1 = curvedsky.rec_lens.qtb(p.lmax,p.rlmin,p.rlmax,lcl[3,:],alm1,blm2,gtype='k',nside=p.nsidet)
      glm2, clm2 = curvedsky.rec_lens.qtb(p.lmax,p.rlmin,p.rlmax,lcl[3,:],alm2,blm1,gtype='k',nside=p.nsidet)

    if q=='EE':
      glm1, clm1 = curvedsky.rec_lens.qee(p.lmax,p.rlmin,p.rlmax,lcl[1,:],alm1,alm2,gtype='k',nside=p.nsidet)
      glm2, clm2 = curvedsky.rec_lens.qee(p.lmax,p.rlmin,p.rlmax,lcl[1,:],alm2,alm1,gtype='k',nside=p.nsidet)

    if q=='EB':
      glm1, clm1 = curvedsky.rec_lens.qeb(p.lmax,p.rlmin,p.rlmax,lcl[1,:],alm1,blm2,gtype='k',nside=p.nsidet)
      glm2, clm2 = curvedsky.rec_lens.qeb(p.lmax,p.rlmin,p.rlmax,lcl[1,:],alm2,blm1,gtype='k',nside=p.nsidet)

    return glm1+glm2, clm1+clm2

  if p.qtype=='rot':
    if q=='EB':
      rlm1 = curvedsky.rec_rot.qeb(p.lmax,p.rlmin,p.rlmax,lcl[1,:],alm1,blm2,nside=p.nsidet)
      rlm2 = curvedsky.rec_rot.qeb(p.lmax,p.rlmin,p.rlmax,lcl[1,:],alm2,blm1,nside=p.nsidet)

    return rlm1+rlm2, (rlm1+rlm2)*0.


@profile
def mean(p,fquad,r):

  for q in p.qlist:

    print('load data first',q)
    glm = np.zeros((p.snmf,p.lmax+1,p.lmax+1),dtype=np.complex)
    clm = np.zeros((p.snmf,p.lmax+1,p.lmax+1),dtype=np.complex)

    for I in range(1,p.snmf+1):
      print(I)
      glm[I-1,:,:], clm[I-1,:,:] = pickle.load(open(fquad[q].alm[I],"rb"))

    print('compute mean field')

    for i in range(p.snmin,p.snmax):
      print(i)
      mfg = np.average(glm,axis=0)
      mfc = np.average(clm,axis=0)
      if i!=0: 
        mfg -= glm[i-1,:,:]/p.snmf
        mfc -= clm[i-1,:,:]/p.snmf
        mfg *= p.snmf/(p.snmf-1.)
        mfc *= p.snmf/(p.snmf-1.)

      print('save to file')
      pickle.dump((mfg,mfc),open(fquad[q].mfb[i],"wb"),protocol=pickle.HIGHEST_PROTOCOL)

      # compute mf cls
      print('cl')
      cl = np.zeros((2,p.lmax+1))
      cl[0,:] = curvedsky.utils.alm2cl(p.lmax,mfg)/r.w4
      cl[1,:] = curvedsky.utils.alm2cl(p.lmax,mfc)/r.w4
      np.savetxt(fquad[q].ml[i],np.concatenate((r.eL[None,:],cl)).T)


