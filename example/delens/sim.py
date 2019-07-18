#* Linear template delensing

# load modules
import numpy as np
import basic
import pickle
import curvedsky

# define parameters
Tcmb = 2.726e6    # CMB temperature
lmax = 3000       # maximum multipole of output cl
rlmin = 500
rlmax = 2048      # reconstruction multipole range
dlmin = 2
dlmax = 2048 
nside = 4096
npix = 12*nside**2
mcnum = 1

# load unlensed and lensed Cls
ucl = basic.aps.read_cambcls('../data/unlensedcls.dat',2,lmax,5)/Tcmb**2
lcl = basic.aps.read_cambcls('../data/lensedcls.dat',2,lmax,4,bb=True)/Tcmb**2

cls = np.zeros((mcnum,7,lmax+1))

for i in range(mcnum):

    print(i)

    # generate gaussian phi
    try:
        plm = pickle.load(open('tmp/phi'+str(i)+'.pkl',"rb"))
    except:
        plm = curvedsky.utils.gauss1alm(lmax,ucl[3,:])
        pickle.dump((plm),open('tmp/phi'+str(i)+'.pkl',"wb"),protocol=pickle.HIGHEST_PROTOCOL)


    # lensed CMB alms
    try:
        Trlm, Erlm, Brlm = pickle.load(open('tmp/lcmb'+str(i)+'.pkl',"rb"))
    except:
        Talm, Ealm = curvedsky.utils.gauss2alm(lmax,ucl[0,:],ucl[1,:],ucl[2,:])
        beta = curvedsky.delens.phi2grad(npix,lmax,plm)
        Trlm, Erlm, Brlm = curvedsky.delens.remap_tp(npix,lmax,beta,np.array((Talm,Ealm,0*Ealm)))
        pickle.dump((Trlm,Erlm,Brlm),open('tmp/lcmb'+str(i)+'.pkl',"wb"),protocol=pickle.HIGHEST_PROTOCOL)

    # aps
    cls[i,0,:] = curvedsky.utils.alm2cl(lmax,Erlm)
    cls[i,1,:] = curvedsky.utils.alm2cl(lmax,Brlm)
    cls[i,2,:] = curvedsky.utils.alm2cl(lmax,plm)

    # lens norm
    try:
        Ag0, Ag1 = np.loadtxt('tmp/al.dat',unpack=True,usecols=(1,2))
    except:
        Ag0, Ac = curvedsky.norm_lens.qtt(lmax,rlmin,rlmax,lcl[0,:rlmax+1],lcl[0,:rlmax+1])
        Ag1, Ac = curvedsky.norm_lens.qeb(lmax,rlmin,rlmax,lcl[1,:rlmax+1],lcl[1,:rlmax+1],lcl[2,:rlmax+1])
        np.savetxt('tmp/al.dat',np.array((np.linspace(0,lmax,lmax+1),Ag0,Ag1)).T)

    # reconstruction
    try:
        plm0, plm1 = pickle.load(open('tmp/qrec'+str(i)+'.pkl',"rb"))
    except:
        # simple diagonal c-inverse
        Fl = np.zeros((3,rlmax+1))
        for l in range(rlmin,rlmax):
            Fl[:,l] = 1./lcl[:3,l]
        # diagonal filtering (since idealistic)
        fTlm = Trlm[0:rlmax+1,0:rlmax+1]*Fl[0,:,None]
        fElm = Erlm[0:rlmax+1,0:rlmax+1]*Fl[1,:,None]
        fBlm = Brlm[0:rlmax+1,0:rlmax+1]*Fl[2,:,None]
        plm0, clm = curvedsky.rec_lens.qtt(lmax,rlmin,rlmax,lcl[0,:rlmax+1],fTlm,fTlm)
        plm1, clm = curvedsky.rec_lens.qeb(lmax,rlmin,rlmax,lcl[1,:rlmax+1],fElm,fBlm)
        plm0 *= Ag0[:,None] 
        plm1 *= Ag1[:,None] 
        pickle.dump((plm0,plm1),open('tmp/qrec'+str(i)+'.pkl',"wb"),protocol=pickle.HIGHEST_PROTOCOL)

    cls[i,3,:] = curvedsky.utils.alm2cl(lmax,plm0)
    cls[i,4,:] = curvedsky.utils.alm2cl(lmax,plm1)

    # template lensing B-mode
    Wl = np.zeros((2,lmax+1))
    for l in range(dlmin,dlmax+1):
        Wl[0,l] = ucl[3,l]/(ucl[3,l]+Ag0[l])
        Wl[1,l] = ucl[3,l]/(ucl[3,l]+Ag1[l])

    plm0 *= Wl[0,:,None]
    plm1 *= Wl[1,:,None]
    blm0 = curvedsky.delens.lensingb(lmax,dlmin,dlmax,dlmin,dlmax,Erlm[:dlmax+1,:dlmax+1],plm0[:dlmax+1,:dlmax+1])
    blm1 = curvedsky.delens.lensingb(lmax,dlmin,dlmax,dlmin,dlmax,Erlm[:dlmax+1,:dlmax+1],plm1[:dlmax+1,:dlmax+1])
    #blm0 = curvedsky.delens.lensingb(lmax,dlmin,dlmax,dlmin,dlmax,Erlm[:dlmax+1,:dlmax+1],plm[:dlmax+1,:dlmax+1],nside=3000)
    #blm1 = blm0

    # aps
    cls[i,5,:] = curvedsky.utils.alm2cl(lmax,blm0)
    cls[i,6,:] = curvedsky.utils.alm2cl(lmax,blm1)

L = np.linspace(0,lmax,lmax+1)
np.savetxt('tmp/resbb.dat',np.concatenate((L[None,:],np.mean(cls,axis=0))).T)

