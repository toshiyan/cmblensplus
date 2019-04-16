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
function rotate_coord, in_uvec, inco=in_coord, outco=out_coord, euler_matrix=eul_mat, stokes_parameters = stokes, help=help
;+
; NAME:
;   rotate_coord
;
;
; PURPOSE:
;   apply change of astronomical coordinates OR active rotation to
;   position vectors and Stokes parameters Q and U
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;    Result= rotate_coord(Vec, Euler_Matrix=, Help=, Inco=, Outco=, Stokes_parameters=)
;
; 
; INPUTS:
;    Vec : a (N x 3) array of coordinates
;
;
; OPTIONAL INPUTS:
;    Help:  if set, prints the documentation header and exits
;
;    Euler_matrix : 3x3 matrix of the active rotation
;      (see eg, Euler_matrix_new routine)
;
;    Inco : input coordinate system (must be used with Outco)
;     'G' :     Galactic
;     'E' :     Ecliptic
;     'C','Q' : Celestial, aka eQuatorial
;
;    Outco : output coordinate system (must be used with Inco)
;      one of 'G', 'E', 'C', 'Q'
;
;    Stokes_parameters : (N x 2) array, values of Stokes parameters Q
;     and U at the location of the input vector vec.
;     where Q measures the polarisation along or across the meridian, 
;     and U the polarisation at 45deg of the meridian
;     On output contain the value of Q and U resulting from the rotation
;      
;
; OUTPUTS:
;    Result : a (N x 3) array of coordinates
;
;
; OPTIONAL OUTPUTS:
;   Stokes_parameters is modified on output
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;    Only one of the rotations (astrophysical coordinate system or
;    active rotation) can be applied at a time
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;    ???: EH, creation
;    2014-08-25: added Help keyword
;-

routine = 'rotate_coord'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if (n_params() ne 1) then begin
    message,'Invalid number of arguments',/noprefix,/infor
    message,'result= '+routine+'(Vec, [Euler_Matrix=, HELP=, Inco=, Outco=, Stokes_parameters=])',/noprefix,/noname

endif

if ((defined(in_coord) or defined(out_coord)) and defined(eul_mat)) then begin
    message,' can not apply the active rotation and the change of coordinate system at the same time',/noprefix
endif

vec0 = transpose([0.d0,0.d0,1.d0]) ; north pole of initial coordinate system

nv = n_elements(in_uvec[*,0])
if (defined(in_coord) and defined(out_coord)) then begin
    ic = strupcase(strmid(in_coord, 0,1)) ; 1st letter, uppercase
    oc = strupcase(strmid(out_coord,0,1))

    loadsky
    if ic eq 'C' then ic = 'Q'
    if oc eq 'C' then oc = 'Q'
    vec = skyconv(in_uvec, inco=ic, outco=oc)

    vec0 = skyconv(vec0, inco=oc, outco=ic)

endif else begin
    vec = in_uvec

    vec0 = vec0

endelse

if (defined(eul_mat)) then begin
    sz = size(eul_mat)
    if ((sz[1] ne sz[2]) or (sz[0] ne 2) or (sz[1] ne 3)) then begin
        message,'The Euler matrix should be 3x3, it is : '+string(sz[0:2],form='(3(i5))'),/noprefix
    endif

    vec  = vec  # eul_mat
    vec0 = vec0 # transpose(eul_mat)

endif

if (defined(stokes)) then begin
    ns = n_elements(stokes[*,0])
    if (ns ne nv) then begin
        message,'Number of Stokes parameter pairs (Q,U) should match number of 3D vectors',/noprefix,/inform
        message,string([ns,nv],form='(2(i10))'),/noprefix
    endif
    
    ; p = initial North Pole
    ; e_phi   = p X r       / sqrt(1-z^2) : unit vector along parallel
    ; e_theta = (p X r) X r / sqrt(1-z^2) : unit vector along meridian
    ; R(e_phi), R(e_theta) : vectors after rotation
    ; R(e) . p = e . Rinv(p)
    ;
    ; R(r) = (x',y',z')
    ; Rinv(p) = (x0', y0', z0')
    ;  abs( p X r ) = sqrt(1-z*z)
    ;   Rinv(p) . r = p . R(r) = z'
    ;
    ; sin(psi) = (p X r) . Rinv(p)      / sqrt(1-z^2) * C
    ;          = -y*x0' + x*y0'         / sqrt(1-z^2) * C
    ; cos(psi) = ((p X r) X r). Rinv(p) / sqrt(1-z^2) * C
    ;          = (z'z - z0')            / sqrt(1-z^2) * C

    ; not normalised SIN(psi) and COS(psi)
    sin_psi = vec0[0,1]*in_uvec[*,0] - vec0[0,0]*in_uvec[*,1]
    cos_psi = - vec0[0,2] + vec[*,2]*in_uvec[*,2]

    ; compute SIN(2*psi) and COS(2*psi)
;     --- slower method, more accurate ?, less memory
;     twopsi = 2.*atan(sin_psi, cos_psi)
;     sin_psi = 0. & cos_psi = 0. ; free memory
;     s2psi = sin(twopsi)
;     c2psi = cos(twopsi)
;     twopsi = 0. ; free memory
;   --- faster method, less accurate ?, more memory
    norm = sin_psi * sin_psi + cos_psi * cos_psi
    s2psi =         2.0d0 * sin_psi * cos_psi / norm
    cos_psi = 0.  ; free memory
    c2psi = 1.0d0 - 2.0d0 * sin_psi * sin_psi / norm
    sin_psi = 0.  ; free memory
    norm    = 0.  ; free memory

    ; new Q and U
    q = stokes[*,0] * c2psi + stokes[*,1] * s2psi
    u = stokes[*,1] * c2psi - stokes[*,0] * s2psi

    stokes = [ [q], [u] ]
    q = 0. & u = 0. ; free memory

endif


return,vec
end

