;+
; NAME:
;       HPX_LATEXIFY
;
;PURPOSE
;	to use PSfrag to do latex syntax correctly
;SYNTAX
;	hpx_latexify, filename, tags, tex, [scale, HEIGHT=, HELP=, OUTNAME=, WIDTH=]
;INPUTS
;	filename: name of eps file
;	tags: name of dummy placeholder
;	tex: the tex you want to replace it
;	scale: the scale you want the latex to be [default to 1]
;NOTES:
;	requires the following latex packages properly installed:
;		geometry, graphicx, psfrag, color, helvet
;	requires you to use font=0
;
;	Follows from discussion by Saurav
;		http://www.astrobetter.com/idl-psfrag/
;	Writes and deletes files called $IDL_TMPDIR/psfrag_*_temp*
;EXAMPLE
;
;   SET_PLOT, 'PS'
;   !P.FONT=0
;   DEVICE, FILENAME='figure.eps'
;   plot,findgen(10),xtitle='xtitle'
;   DEVICE, /CLOSE
;
;   hpx_latexify, 'figure.eps', 'xtitle', '$M_{\odot}$'
;	
; HISTORY:
; adapted from latexify.pro by R. da Silva, UCSC, 1-4-11
;  and from LS_latexify.pro by L. Spencer (2011)
;  version 1.0: 2015-05-12
;	
;-
pro hpx_latexify, filename, tag, tex, scale, $
                  height = height, $
                  help = help, $
                  outname=outname, $
                  width = width

routine = 'hpx_latexify'
syntax = routine+', filename, tags, tex, [scale, HEIGHT=, HELP=, OUTNAME=, WIDTH=]'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if n_params() lt 3 then begin
    print,syntax
    return
endif

if ~keyword_set(scale) then scale=replicate(1., n_elements(tag))
scale=strtrim(scale*1.,2)

; temporary files
tmpdir = getenv('IDL_TMPDIR')
stime  = string(long(systime(1)*1e3 mod 1.e8), form='(i8.8)')
tmprad = add_prefix(tmpdir,path_sep(),/suff)+'psfrag_'+stime+'_temp'
noname=0
if ~keyword_set(outname) then begin
    outname=tmprad + '2.ps'
    noname=1
endif

; generate LaTeX file
openw, lun, tmprad + '.tex', /get_lun
printf, lun, '\documentclass{article}'
printf, lun, '\nonstopmode' ; will crash on error instead of prompting user
printf, lun, '\usepackage{geometry, graphicx, psfrag}'
printf, lun, '\usepackage{color}'
; printf, lun, '\usepackage{helvet}'
; printf, lun, '\renewcommand{\familydefault}{\sfdefault}'
printf, lun,'\pagestyle{empty}'
printf, lun,'\geometry{paperwidth='+strtrim(width*1.01,2)+'cm,'+$
                    'paperheight='+strtrim(height*1.01,2)+'cm,margin=0pt}'
printf, lun,'\begin{document}'
for i=0, n_elements(tag)-1 do begin
    printf, lun,'\psfrag{'+tag[i]+'}[cc][cc]['+scale[i]+']{'+tex[i]+'}'
endfor
printf, lun,'\begin{center}'
printf, lun,'\includegraphics[width='+$
    	strtrim(width,2)+'cm,height='+strtrim(height,2)+'cm]{'+filename+'}'
printf, lun,'\end{center}'
printf, lun,'\end{document}'
close, lun
free_lun, lun

; generate PS file
spawn, 'cd '+tmpdir+'; latex '               +tmprad +'.tex', ls0
if ~file_test(tmprad+'.dvi') then begin
    print, ls0,form='(a)'
    message,'Error in LaTeX. Aborting.'
endif else begin
    command = 'cd '+tmpdir+'; dvips -o '+outname+' '+tmprad +'.dvi'
    if (strupcase(!version.os_family) eq 'WINDOWS') then begin
        spawn,      command, ls1
    endif else begin
        spawn, /sh, command+' > '+tmprad+'_dvi.log  2>&1 '
    endelse
; spawn, 'dvips -o '+tmprad+'.ps '+tmprad+'.dvi', ls1
; spawn, 'ps2epsi  '+tmprad+'.ps '+tmprad+'.epsi, ls2
; spawn, "perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/' <"+$
;	 tmprad + '.epsi > '+outname, ls3
endelse

;garbage cleanup
;;;;;;; spawn,'cp '+filename+' '+tmprad+'_in.ps'
if noname then begin
	spawn, 'mv -f '+outname+' '+filename
	outname=filename
endif

;;;;;;; print,tmprad+['2.ps','.tex','.aux','.log','.dvi','_dvi.log']
;;;print,'Files NOT deleted: ',tmprad+['2.ps','.tex','.aux','.log','.dvi','_dvi.log']
file_delete, tmprad+['2.ps','.tex','.aux','.log','.dvi','_dvi.log'], allow_nonexistent=1


end
