; sample initialization file for stokes3
; SC lattice config of 8 particles in (5,5,5) box
; $Id: stokes3.scm,v 1.1 2006/10/23 17:15:09 kichiki Exp $

(define outfile    "stokes3.SC8.nc") ; output filename

(define version    "F")     ; version. "F", "FT", or "FTS"
(define flag-mat   #t)      ; #t => matrix scheme, #f => atimes scheme
(define flag-lub   #t)      ; #t => with lub,      #f => without lub

(define np         8)       ; number of particles
(define nm         8)       ; number of free particles
(define nloop      1000)    ; main loop
(define dt         1.0e-1)  ; time interval
(define stokes     0.0)     ; effective stokes number
(define ewald-tr   4.1)     ; time ratio Tr/Tk for Ewald summation
(define ewald-eps  1.0e-12) ; cut-off limit for Ewald summation

(define Ui '(0.0  0.0  0.0)) ; imposed translational velocity
(define Oi '(0.0  0.0  0.0)) ; imposed angular velocity relative to O
(define Ei '(0.0  0.0  0.0  0.0  0.0)) ; imposed strain relative to O

(define F0 '(0.0  0.0 -0.1)) ; applied force
(define T0 '(0.0  0.0  0.0)) ; applied torque

; lattice vector
(define lattice '(5.0  5.0  5.0))

; initial configuration
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
