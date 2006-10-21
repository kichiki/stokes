; sample set file for xi3
; SC lattice config of 8 particles in (5,5,5) box
; $Id: xi3.scm,v 1.1 2006/10/21 19:05:27 kichiki Exp $
(define version    "F")     ; version. "F", "FT", or "FTS"
(define flag-mat   t)       ; t => matrix scheme, nil => atimes scheme
(define flag-notbl nil)     ; t => no-table,      nil => with table

(define np         8)       ; number of particles
(define ewald-eps  1.0e-12) ; cut-off limit for Ewald summation

; lattice vector
(define lattice '(5.0 5.0 5.0))

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

; list of time ratio Tr/Tk for Ewald summation (optional)
;(define ewald-trs
;  '(0.1
;    1.0
;    10.0
;    100.0
;    ))
