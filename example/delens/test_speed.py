# * Linear template delensing

# load modules
import numpy as np
import basic
import curvedsky
import time

# parameters
lmin = 1        # minimum multipole
lmax = 2000     # maximum multipole

# input phi and E
iplm = np.load('/project/projectdirs/sobs/delensing/b_template_tests/input_phi_alm_lmax2000.npy')
iElm = np.load('/project/projectdirs/sobs/delensing/b_template_tests/e_len_alm_lmax2000.npy')

# convert healpy to healpix alm
print(len(iplm))
Elm = curvedsky.utils.lm_healpy2healpix(len(iElm),iElm,lmax)
plm = curvedsky.utils.lm_healpy2healpix(len(iplm),iplm,lmax)

# speed of template lensing B-mode
start_time = time.time()
lalm = curvedsky.delens.lensingb(lmax,lmin,lmax,lmin,lmax,Elm,plm)
end_time = time.time()

print('time taken: '+str(end_time-start_time)+' s')


