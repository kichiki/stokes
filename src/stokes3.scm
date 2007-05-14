; sample initialization file for stokes3
; SC lattice config of 8 particles in (5,5,5) box
; $Id: stokes3.scm,v 1.5 2007/05/14 07:57:03 kichiki Exp $

;; output parameters
(define outfile    "stokes3.SC8.nc") ; output filename
(define dt         1.0)     ; time interval
(define nloop      1000)    ; main loop
(define flag-Q     #f)      ; #t => output quaternion, #f => no quaternion.

;; core libstokes parameters
(define version    "F")     ; version. "F", "FT", or "FTS"
(define flag-mat   #t)      ; #t => matrix scheme, #f => atimes scheme
(define flag-lub   #t)      ; #t => with lub,      #f => without lub
(define lub-min    2.0000000001) ; min cutoff of distance for lub
(define lub-max    4.0)          ; max cutoff of distance for lub

;; periodic systems
(define periodic   #t)      ; #f => non periodic, #t => periodic
(define ewald-tr   4.1)     ; time ratio Tr/Tk for Ewald summation
(define ewald-eps  1.0e-12) ; cut-off limit for Ewald summation
(define lattice    '(5.0  5.0  5.0)) ; size of the periodic box

;; ODE parameters
(define ode-solver "rkf45")
; the following solvers are available
; "rk2"       Embedded Runge-Kutta (2, 3) method.
; "rk4"       4th order (classical) Runge-Kutta.
; "rkf45"     Embedded Runge-Kutta-Fehlberg (4, 5) method.
; "rkck"      Embedded Runge-Kutta Cash-Karp (4, 5) method.
; "rk8pd"     Embedded Runge-Kutta Prince-Dormand (8,9) method.
; "rk2imp"    Implicit 2nd order Runge-Kutta at Gaussian points.
; "rk4imp"    Implicit 4th order Runge-Kutta at Gaussian points.
; "gear1"     M=1 implicit Gear method.
; "gear2"     M=2 implicit Gear method.
(define ode-eps 1.0e-6) ; ODE control parameter eps

;; system parameters
(define np         8)       ; number of particles
(define nm         8)       ; number of free particles

; particle radius
; set '() for monodisperse (a = 1.0 for all particles)
;(define a '())
; otherwise, poly codes are used in the calculations
(define a #(
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
))

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

(define Ui '(0.0  0.0  0.0)) ; imposed translational velocity
(define Oi '(0.0  0.0  0.0)) ; imposed angular velocity relative to O
(define Ei '(0.0  0.0  0.0  0.0  0.0)) ; imposed strain relative to O

(define F0 '(0.0  0.0 -0.1)) ; applied force
(define T0 '(0.0  0.0  0.0)) ; applied torque

(define stokes 0.0)     ; effective stokes number
(define ncol   10)      ; frequency of collision check in dt for st != 0


;; bond parameters
(define bonds '())
(define flag-relax #f) ; #f => stokesian dynamics with bond interactions
                       ; #t => relaxation dynamics only with bond interactions
(define gamma 1.0)     ; friction coefficient for relaxation dynamics
