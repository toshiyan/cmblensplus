#! /usr/bin/env python3

from __future__ import print_function # so that print outputs look the same in python 2 and 3
#
import sys, getopt, os, time
import healpy          as hp
from datetime  import datetime

#
try:
    import astropy.io.fits as pf
except:
    import pyfits as pf
#
try:
    import pytz
    UTC = pytz.timezone('UTC')
except:
    from astropy.time import TimezoneInfo
    UTC = TimezoneInfo()

'''
  change_polcconv
   routine to change the polarization coordinate convention in
   HEALPix map FITS files.

   Changing this convention will change the sign of the U Stokes parameter,
   and swap the POLCCONV FITS keyword between "IAU" and "COSMO"

   Can be called from within python (2 and 3):
      change_polcconv(file_in, file_out, ctype, force=False)

   or as a shell command:
      > python change_polcconv.py [--force] [--c2c|--c2i|--i2c|--i2i] file_in file_out
   or simply (if python3 is available):
      > ./change_polcconv.py [--force] [--c2c|--c2i|--i2c|--i2i] file_in file_out

   with:
      file_in:  input  HEALPix map FITS file

      file_out: output HEALPix map FITS file, must be different from file_in

      force: if set, the value of POLCCONV read from input file is ignored.
           The sign of U is swapped (if used with c2i or i2c),
           and/or the keyword is updated to IAU (if used with i2i and c2i)
           or to COSMO (if used with c2c or i2c).

      ctype must be among "c2c", "c2i", "i2c", "i2i" (NOT case sensitive)

          c2c: does NOT change coordinate system
              -if POLCCONV is found with value "IAU", and force is not set,
               program will issue error message and no file is written;
              -in all other cases POLCCONV is set/reset to "COSMO", but
               data is NOT changed

          c2i: changes from COSMO to IAU coordinate convention
              -if POLCCONV is not found or found with value "COSMO",
               or force is set, it is added/replaced with value "IAU",
               and the sign of the U stokes parameter map is changed;
              -if POLCCONV already has value "IAU",
               and force is NOT set, 
               file_in is copied unchanged into file_out

          i2c: changes from IAU to COSMO coordinate convention
              -if POLCCONV is not found or found with value "IAU",
               or force is set, it is added/replaced with value "COSMO",
               and the sign of the U stokes parameter map is changed;
              -if POLCCONV already has value "COSMO",
               and force is NOT set, 
               file_in is copied unchanged into file_out

          i2i: does NOT change coordinate system
              -if POLCCONV is found with value "COSMO", and force is not set,
               program will issue error message and no file is written;
              -in all other cases POLCCONV is set/reset to "IAU", but
               data is NOT changed


   HISTORY: 2018-05: v1.0, E. Hivon
                     adapted from HEALPix/IDL change_polcconv.pro
            2020-02: v1.1, process partial P,T,Q,U files

'''

#--------------------------------

def fitstype(hdr_or_file, ext=1):
    '''ftype = fitstype(Header or FITS file, ext=1)
       returns Healpix FITS type
       0: no extension
       1: ASCII extension (mostly for C(l) power spectra)
       2: binary extension with implicit pixel indexing (for full sky data)
       3: binary extension with explicit pixel indexing (for cut sky data)
       '''
    
    if isinstance(hdr_or_file, pf.Header):
        header = hdr_or_file
    else:
        header = pf.getheader(hdr_or_file, ext=ext)

    undef = 'UNDEF'

    xt = (header.get('XTENSION', default=undef)).upper()
    if xt.startswith(undef):     ft = 0
    if xt.startswith('TABLE'):   ft = 1
    if xt.startswith('BINTABLE'):ft = 2

    # adapted from healpy 1.11.0   read_map
    #
    # partial sky: check OBJECT, then INDXSCHM
    partial = False
    if (ft == 2):
        obj = (header.get('OBJECT', default=undef)).strip()
        if obj != undef:
            if obj == 'PARTIAL':
                partial = True
            elif obj == 'FULLSKY':
                partial = False
        
        schm = (header.get('INDXSCHM', default=undef)).strip()
        if schm != undef:
            if schm == 'EXPLICIT':
                if obj == 'FULLSKY':
                    raise ValueError('Incompatible INDXSCHM keyword')
                partial = True
            elif schm == 'IMPLICIT':
                if obj == 'PARTIAL':
                    raise ValueError('Incompatible INDXSCHM keyword')
                partial = False
        
        if schm == undef:
            schm = (partial and 'EXPLICIT' or 'IMPLICIT')
            # print("No INDXSCHM keyword in FITS header, assume %s"%(schm))
        
        if partial: ft = 3
    
    return ft

#--------------------------------

def change_polcconv(file_in, file_out, ctype, force=False):

    t0 = time.time()
    if not os.path.exists(file_in):
        sys.exit('Input file %s does not exist'%(file_in))

    if (file_in == file_out):
        #raise ValueError('Input and output files must be different')
        sys.exit('Input and output files must be different')
    
    uctype = (ctype.upper())[2:5]
    allowed = ('C2C','C2I','I2C','I2I')
    if (not uctype.startswith(allowed)):
        print ('ERROR: Ctype must be among: ')
        print (allowed)
        print ('While it is: %s'%(uctype))
        #raise ValueError('Aborting')
        sys.exit('Aborting')
    
    do_c2c =uctype.startswith('C2C')
    do_c2i =uctype.startswith('C2I')
    do_i2c =uctype.startswith('I2C')
    do_i2i =uctype.startswith('I2I')

    kw_pol = 'POLCCONV'
    kw_iau = 'IAU'
    kw_cos = 'COSMO'
    undef  = 'UNDEF'
    modify_map = (do_c2i or do_i2c)



    print ('%s -> [%s] -> %s'%(file_in, uctype, file_out))
    # file type
    info = pf.info(file_in, output=False)
    hdr  = pf.getheader(file_in, ext=1)
    ftype = fitstype(hdr)
    n_ext = len(info)-1

    # copy primary header
    hdr0      = pf.getheader(file_in, ext=0)
    polcconv1 = hdr0.get(kw_pol, default=undef).upper()
    fullout   = pf.HDUList([pf.PrimaryHDU(data=None, header=hdr0)])
    
    # loop on extensions
    for ext in range(1,n_ext+1):
        extname = info[ext][1]
        if (ftype == 3): # cut sky
            dcut, xhdr = pf.getdata(file_in, ext=ext, header=True)
            names = dcut.names
        else: # full sky
            tqu, xhdr  = pf.getdata(file_in, ext=ext, header=True)
            names = tqu.names

        nside = xhdr.get('NSIDE')
        in_xxx = True ; in_iau = False ; in_cos = False
        str_ext = 'extension #%s (%s)'%(str(ext),extname)
        nmaps = len(names)
        xxQU= (nmaps >= 4 and     # first 4 column names are *,*,Q*,U*
               names[2].startswith(('Q','q')) and
               names[3].startswith(('U','u')) )
        xQU = (nmaps >= 3 and     # first 3 column names are *,Q*,U*
               names[1].startswith(('Q','q')) and
               names[2].startswith(('U','u')) )
        QU  = (nmaps == 2 and     # first 2 column names are Q*,U*
               names[0].startswith(('Q','q')) and
               names[1].startswith(('U','u')) )

        # check POLCCONV in input file
        polcconv = xhdr.get(kw_pol, default=undef).upper()
        if (polcconv1 != undef and polcconv != undef and polcconv1 != polcconv):
            raise ValueError('Unconsistent %s values found in Primary and extension headers'%(kw_pol))
        #in_xxx = (polcconv == undef)
        in_iau = (polcconv == kw_iau)
        in_cos = (polcconv == kw_cos)
        in_xxx = not (in_iau or in_cos)
        done    = False
        non_pol = False
        do_edit = False
        head_only = False

        if not force:
            if ((in_cos and do_i2i) or (in_iau and do_c2c)):
                print ('FITS file: %s'%(file_in))
                print ('%s: %s'%(kw_pol, polcconv))
                print ('Can not perform requested %s edition'%(uctype))
                print ('Use the force=True keyword to only update %s value.'%(kw_pol))
                sys.exit(1)

            done = (in_cos and (do_i2c or do_c2c) or in_iau and (do_c2i or do_i2i))

        if (modify_map and not done):
        # change sign of U and related quantities
            clist = [] ; do_edit = False
            if (ftype == 3):
                if (xxQU):           # change sign of U (4th column of P,T,Q,U)
                    clist = [3] ; do_edit = True
                else:
                    if (ext == 3): # change sign of U (2nd column of P,S,N,S; last extension)
                        clist = [1] ; do_edit = True
                    else:
                        head_only = (n_ext == 3) # update header of 1st and 2nd extensions
            else:
                if (nmaps <= 1):
                    non_pol = True

                if (nmaps == 2):
                    if (n_ext == 1 and QU):
                        # Planck 2 columns (Q*, U*)
                        clist = [1] ; do_edit = True
                    else:
                        non_pol = True

                if (nmaps == 3):
                    # standard Healpix full sky format (T, Q, U ; TT, QQ, UU ; QU, TU, TQ)
                    if (ext == 1):
                        if (xQU):
                            clist = [2] ; do_edit = True # change sign of U
                        else:
                            non_pol = True
                    if (ext == 3): # change sign of QU and of TU
                        clist = [0,1] ; do_edit = True

                if (n_ext >= 1 and nmaps == 5):
                    # Planck 5 cols (I_STOKES, Q_STOKES, U_STOKES, TMASK, PMASK) + BEAM_TF extension
                    if (ext == 1 and xQU):
                        clist = [2] ; do_edit = True # change sign of U

                if (n_ext == 1 and nmaps == 10):
                    # Planck R3 (2017) 10 cols (I_STOKES, Q_STOKES, U_STOKES, HIT, II_COV, IQ_COV, IU_COV, QQ_COV, QU_COV, UU_COV)
                    if (ext == 1 and xQU):
                        clist = [2, 6, 8] ; do_edit = True # change sign of U, IU, QU

                if (nside <= 512 and n_ext == 2):
                    # WMAP: https://lambda.gsfc.nasa.gov/product/map/dr5/skymap_file_format_info.cfm
                    if (nmaps >= 5):
                        # WMAP (I,Q,U,S) format (T, Q, U, S, N_OBS ; 
                        #      N_OBS, M11=SS, M12=SQ, M13=SU, M23=QU, M22=QQ, M33=UU)
                        if (ext == 1 and nmaps == 5 and xQU):
                            clist = [2] ; do_edit = True # change sign of U
                        if (ext == 2 and nmaps == 7):
                            clist = [3,4] ; do_edit = True # change sign of SU and QU
                            
                    if (nmaps == 4): # WMAP (I,Q,U) format (T, Q, U, N_OBS ; N_OBS, QQ, QU, UU)
                        clist = [2] ; do_edit = True # change sign of U and of N_QU
                        #                 else:
                        #                     raise ValueError('Unrecognised format in %s'%file_in)
                if (not do_edit):
                    if (non_pol):
                        print (' Unpolarised format, no change in %s'%(str_ext))
                    else:
                        print ('WARNING: unrecognised format in %s'%(str_ext))
                        print (nmaps, n_ext)
                        print (names)
                        done = True
                        
            if (do_edit):
                xhdr.add_history(file_in)
                xhdr.add_history('edited on %s UTC'%(datetime.isoformat(datetime.now(UTC))[0:19]))
                for c in clist:
                    if (ftype == 3):
                        dcut[names[c]] *= -1
                    else:
                        tqu[names[c]] *= -1
                    str_col = 'column #%s (%s)'%(str(c+1),names[c].strip())
                    print (' Flip sign of %s in %s'%(str_col,str_ext))
                    xhdr.add_history('modified %s to match coord. conv. (%s)'%(str_col,kw_pol))
                
        if ( (do_edit or in_xxx or force or head_only) and not non_pol):
            new_cconv = kw_cos if (do_c2c or do_i2c) else kw_iau
            xhdr[kw_pol]= (new_cconv, ' Coord. convention for polarisation (COSMO/IAU)')
            print (' %s keyword %s=%s in header of %s'%((in_xxx and 'Add' or 'Update'),kw_pol,new_cconv,str_ext))
        else:
            print (' No change in header of %s'%(str_ext))

        # append new data
        if (ftype == 3):
            out = pf.BinTableHDU.from_columns(dcut)
        else:
            out = pf.BinTableHDU.from_columns(tqu)
        out.header = xhdr
        fullout.append(out)

    # write new file
    fullout.writeto(file_out, overwrite=True)
    # list input and output file
    os.system('ls -l %s %s'%(file_in, file_out))
    t1 = time.time()
    print ("Time [s] = ",t1-t0)
    print ('     ----------------        ')
#--------------------------------

if (__name__ == "__main__"):
    
#     change_polcconv('/tmp/tqu.fits','/tmp/tqu_x.fits','x')

#     change_polcconv('/tmp/tqu.fits','/tmp/tqu_IAU.fits','C2I')

#     change_polcconv('/tmp/tqu.fits','/tmp/tqu_COS.fits','C2C')

#     change_polcconv('/tmp/tqu.fits','/tmp/tqu_i2c.fits','I2C')

#     change_polcconv('/tmp/tqu.fits','/tmp/tqu_i2i.fits','I2I',force=True)

#     change_polcconv('/tmp/tqu.fits','/tmp/tqu_i2i.fits','I2I')

    routine = "change_polcconv.py"
    usage1  = "usage: python " + routine + " [--force] [--c2c|--c2i|--i2c|--i2i] file_in file_out"
    allowed = ('c2c','c2i','i2c','i2i')

    try:
        opts, args = getopt.getopt(sys.argv[1:], \
                                   '', \
                                   ['force']+list(allowed))
    except getopt.GetoptError:
        print (usage1)
        sys.exit(2)


    if (len(args) != 2 or len(opts) < 1 or len(opts)> 2):
        print (usage1)
        sys.exit(0)

    file_in  = args[0]
    file_out = args[1]

    force = False
    for opt, val in opts:
        if (opt == '--force'):
            force = True
        elif (opt.startswith(allowed,2)):
            ctype = opt
        else:
            print (usage1)
            sys.exit(2)

    print (file_in, file_out, force, ctype)
    change_polcconv(file_in, file_out, ctype, force=force)

