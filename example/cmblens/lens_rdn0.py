# example of RDN0 bias calculation
import numpy as np
import analytic
import curvedsky
import pickle

snn0 = 50
snrd = 100

lmax = 3000
rlmin = 500
rlmax = 3000
eL = np.linspace(0,lmax,lmax+1)

cltt = 
Fl = 1./cltt

falm = ['alm.dat']
frdn0 = ['rdn0.dat']

# normalization
Alg = {}
Alc = {}
Alg['TT'], Alc['TT'] = analytic.rec_lens.qtt(lmax,rlmin,rlmax,cltt,cltt)

# compute N0
n0 = np.zeros((snn0,2,lmax+1))
for i in range(snn0):
  print (2*i, 2*i+1)

  Talm1 = Fl*pickle.load(open(falm[2*i],"rb"))
  Talm2 = Fl*pickle.load(open(falm[2*i+1],"rb"))

  glm1, clm1 = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,cltt,Talm1,Talm2)
  glm2, clm2 = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,cltt,Talm2,Talm1)
  glm1 *= Alg['TT'][:,None]
  glm2 *= Alg['TT'][:,None]
  clm1 *= Alc['TT'][:,None]
  clm2 *= Alc['TT'][:,None]

  n0[i,0,:] = curvedsky.utils.alm2cl(lmax,glm1+glm2)/2.
  n0[i,1,:] = curvedsky.utils.alm2cl(lmax,clm1+clm2)/2.

n0 = np.mean(n0,axis=0)

# compute RDN0
for i in range(1):
  print(i)

  Talm1 = Fl*pickle.load(open(falm[i],"rb"))

  cl = np.zeros((2,lmax+1))
  sn = 0.
  for I in range(snrd):

    if I == i:  continue

    print(I)
    Talm2 = Fl*pickle.load(open(falm[I],"rb"))

    glm1, clm1 = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,cltt,Talm1,Talm2)
    glm2, clm2 = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,cltt,Talm2,Talm1)
    glm1 *= Alg['TT'][:,None]
    glm2 *= Alg['TT'][:,None]
    clm1 *= Alc['TT'][:,None]
    clm2 *= Alc['TT'][:,None]

    sn += 1.
    cl[0,:] += curvedsky.utils.alm2cl(lmax,glm1+glm2)
    cl[1,:] += curvedsky.utils.alm2cl(lmax,clm1+clm2)

  cl = cl/sn - n0
  print ('save RDN0')
  np.savetxt(frdn0[i],np.concatenate((eL[None,:],cl)).T)

