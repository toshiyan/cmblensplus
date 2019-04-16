; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2013  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.sourceforge.net
;
; -----------------------------------------------------------------------------
pro change_polcconv, file_in, file_out, i2c=i2c, c2i=c2i, c2c=c2c, i2i=i2i, force=force
;+
; NAME:
;          change_polcconv
;
; PURPOSE:
;          change the coordinate convention in FITS file containing a polarised sky map.
;          The main effect is to change the sign of the U Stokes parameter,
;          and add/update the POLCCONV FITS header with either
;          COSMO or IAU value
;
; CATEGORY:
;
; CALLING SEQUENCE:
;          CHANGE_POLCCONV, File_In , File_Out, [/I2C, /C2I, /C2C, /I2I, /FORCE]
;
; INPUTS:
;           File_In: Input FITS file
;
;           File_Out: Output FITS file
;
; KEYWORD PARAMETERS:
;          one and only one of the following four should be set
;
;          I2C: changes from IAU to COSMO coordinate convention
;              -if POLCCONV is not found or found with value 'IAU', it is
;               added/replaced with value 'COSMO', and the sign of the U stokes
;               parameter map is changed
;              -if POLCCONV already has value 'COSMO', File_In is copied
;              unchanged into File_Out
;
;          C2I: changes from COSMO to IAU coordinate convention
;              -if POLCCONV is not found or found with value 'COSMO', it is
;               added/replaced with value 'IAU', and the sign of the U stokes
;               parameter map is changed
;              -if POLCCONV already has value 'IAU',  File_In is copied
;              unchanged into File_Out
;
;          C2C: does NOT change coordinate system
;              -if POLCCONV is found with value 'IAU', program will issue error
;              message and no file is written
;              -in all other case POLCCONV is set/reset to 'COSMO', but
;               data is NOT changed
;
;          I2I: does NOT change coordinate system
;              -if POLCCONV is found with value 'COSMO', program will issue error
;              message and no file is written
;              -in all other case POLCCONV is set/reset to 'IAU', but
;               data is NOT changed
;
;          FORCE: if set, the value of POLCCONV read from the FITS header is
;          ignored.
;              The sign of U is swapped (if used with C2I or I2C), and the
;              keyword is updated accordingly
;
;
;
; OUTPUTS:
;          none
;
; OPTIONAL OUTPUTS:
;          none
;
; COMMON BLOCKS:
;          none
;
; SIDE EFFECTS:
;          may write a new FITS file or modify an existing one
;
; RESTRICTIONS:
;          only works with recognised FITS formats
;
; PROCEDURE:
;          read file, modify header, write file
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;          v1.0, May 2005, Eric Hivon
;          v1.1, Oct 2006, Eric Hivon: added /force
;-


syntax = 'CHANGE_POLCCONV, File_In , File_Out [, /I2C, /C2I, /C2C, /I2I, /FORCE]'

if (n_params() ne 2) then begin
    print,syntax
    return
endif

do_i2c = keyword_set(i2c) 
do_c2i = keyword_set(c2i)
do_c2c = keyword_set(c2c)
do_i2i = keyword_set(i2i)
do_force = keyword_set(force)
sum = do_i2c + do_c2i + do_c2c + do_i2i
if (sum ne 1) then begin
    print,syntax
    message,'One and only one keyword out of /I2C, /C2I, /C2C, /I2I should be set'
endif

; full sky TQU format (with extra extension for correlation)
; WMAP format (extra column and extra extension with correlation)
; cut sky format (multi extension ?)

kw_pol = 'POLCCONV'
kw_iau = 'IAU'
kw_cos = 'COSMO'
fits_info, File_In, /silent, n_ext = n_ext
junk = getsize_fits(File_In, type = type)

if (n_ext ne 3 and type eq 3) then begin
    message,'Unpolarised format cut sky format. No changes made',/info
    return
endif

print,file_in, file_out
for i=0,n_ext-1 do begin
    ; read data
    if (type eq 3) then begin
        read_fits_cut4, File_In, pixel, signal, n_obs, serror, hdr=hdr, xhdr=xhdr, ext = i
    endif else begin
        read_tqu, File_In, tqu, hdr=hdr, xhdr = xhdr, ext = i
    endelse
    xhdr_tmp = [hdr,xhdr]
    polcconv = strupcase(sxpar(xhdr_tmp,KW_POL,count=count))
    in_xxx = 1 & in_iau = 0 & in_cos = 0
    modify_map = (do_c2i or do_i2c)
    
    case count of
        0: in_xxx = 1
        1: begin
            in_iau = (strtrim(polcconv) eq KW_IAU)
            in_cos = (strtrim(polcconv) eq KW_COS)
            in_xxx = 1 - (in_iau + in_cos)
        end
        else: begin
            print,'FITS file: ',strtrim(File_In)
            message,'several '+KW_POL+' keyword in header of extension'+string(i)
        end
    endcase
    if (1 - do_force) then begin
        if ((in_cos and do_i2i) or (in_iau and do_c2c)) then begin
            print,'FITS file: ',strtrim(File_In)
            print,KW_POL+': ',polcconv
            message,'can not perform requested edition'
        endif
        if (in_cos *(do_i2c+do_c2c) or in_iau *(do_c2i+do_i2i)) then goto, skip
    endif

    if (modify_map) then begin
        ; change sign of U and related quantities
        do_edit = 0
        nmaps = n_elements(tqu[0,*])
        if (type eq 3) then begin ; cut sky format: change sign of U (last extension)
            if (i eq 2) then begin
                signal = -signal
                do_edit = 1
            endif
        endif else begin
            case nmaps of 
                4: begin ; WMAP format
                    tqu[*,2] = - tqu[*,2] ; change sign of U and of N_QU
                    do_edit = 1
                end
                3: begin; standard Healpix full sky format
                    if (i eq 0) then begin
                        tqu[*,2] = - tqu[*,2] ; change sign of U
                        do_edit = 1
                    endif
                    if (i eq 2) then begin
                        tqu[*,0:1] = - tqu[*,0:1] ; change sign of QU and of TU
                        do_edit = 1
                    endif
                end
                1 or 2: begin
                    print,'Warning: Unpolarised format, no change in extension '+string(i)
                    goto, skip
                end
                else: begin
                    message,'unrecognised format'
                end
            endcase
        endelse
        if (do_edit) then begin
            sxaddpar, xhdr, 'HISTORY','modified U related components to match coordinate convention (POLCCONV)'
            print,'flip sign of U map in extension '+string(i)
        endif
    endif
    ; write new data
    new_cconv = (do_c2c+do_i2c) ? KW_COS : KW_IAU
    sxaddpar, xhdr, KW_POL, new_cconv,' Coord. convention for polarisation (COSMO/IAU)'
;    sxaddpar,info_header,'COMMENT',' either IAU or COSMO',after=kw_pol
    print,'Update/add '+strtrim(KW_POL)+'in extension '+string(i)
    skip:
    if (type eq 3) then begin
        write_fits_cut4, File_Out, pixel, signal, n_obs, serror, hdr=hdr, xhdr=xhdr, /pol, ext=i
    endif else begin
        write_tqu, File_Out, tqu, hdr=hdr, xhdr=xhdr, ext = i
    endelse
endfor

return
end
