function rodrigues_matrix, vector, angle, deg=deg, help=help
;+
;  Rot_matrix = Rodrigues_matrix(Vector, Angle [, DEG=, HELP=])
;  
;   returns the 3x3 matrix describing the rotation of
;   Angle (in Radians or in Degrees) around Vector, using
;   the Rodrigues formula
;   http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
;
;
; INPUTS:
;
;     Vector: 3-element vector describing the axis of rotation
;         only its normalized version will be used
;
;     Angle: scalar, angle of rotation in Radians,
;            unless DEG is set
;
;
; KEYWORD PARAMETERS:
;
;       DEG : if set, the angle is measured in degrees
;
;       HELP: if set, the documentation header is printed
;
;
; RELATED ROUTINES:
;
;       Euler_matrix_new
;       Rotate_coord
;
; MODIFICATION HISTORY:
;    2014-08-26: EH, creation
;
;
;-

routine='rodrigues_matrix'
syntax = 'rot = '+routine+'(Vector, Angle [, DEG=, HELP=]) '

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() ne 2 then begin
    print,syntax
    return,-1
endif

; normalize vector
v = vector[0:2]
v /= sqrt(total(v*v))

; cross product matrix
kmat = [[ 0,    -v[2],  v[1]],$
        [ v[2],  0,    -v[0]],$
        [-v[1],  v[0],  0   ] ]

; trigonometry
double = SIZE(angle,/TNAME) EQ 'DOUBLE'
convert = 1.0
IF KEYWORD_SET(deg) THEN convert = double ? !dpi/180d : !DtoR
c1 = COS(angle[0]*convert)
s1 = SIN(angle[0]*convert)

; Rotation matrix: Rodrigues formula
rot = diag_matrix([1,1,1]) + s1 * kmat + (1-c1) * kmat#kmat

return, rot
end
