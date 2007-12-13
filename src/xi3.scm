; sample set file for xi3
; SC lattice config of 8 particles in (5,5,5) box
; $Id: xi3.scm,v 1.3 2007/12/13 06:11:28 kichiki Exp $
(define version    "F")     ; version. "F", "FT", or "FTS"
(define flag-mat   #t)      ; #t => matrix scheme, #f => atimes scheme
(define flag-notbl #f)      ; #t => no-table,      #f => with table

(define np         8)       ; number of particles
(define ewald-eps  1.0e-12) ; cut-off limit for Ewald summation

; lattice vector
(define lattice '(5.0  5.0  5.0))

; configuration of particles
(define x #(
0.0  0.0  0.0
2.5  0.0  0.0
0.0  2.5  0.0
0.0  0.0  2.5
0.0  2.5  2.5
2.5  0.0  2.5
2.5  2.5  0.0
2.5  2.5  2.5
))

; particle radius
(define a '())

; slip length
(define slip '())

; list of time ratio Tr/Tk for Ewald summation (optional)
;(define ewald-trs
;  '(0.1
;    1.0
;    10.0
;    100.0
;    ))
