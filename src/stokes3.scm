; sample initialization file for stokes3
; SC lattice config of 8 particles in (5,5,5) box
; $Id: stokes3.scm,v 1.9 2007/12/26 06:45:49 kichiki Exp $

;; output parameters
(define outfile    "stokes3.SC8.nc") ; output filename
(define dt         1.0)     ; time interval
(define nloop      1000)    ; main loop
(define flag-Q     #f)      ; #t => output quaternion, #f => no quaternion.

;; core libstokes parameters
(define version    "F")     ; version. "F", "FT", or "FTS"
(define flag-mat   #t)      ; #t => matrix scheme, #f => atimes scheme
(define flag-lub   #t)      ; #t => with lub,      #f => without lub
(define rmin       0.0)     ; param for min distance (ai+aj) * rmin
(define lub-min    2.0000000001) ; min cutoff of distance for lub
(define lub-max    4.0)          ; max cutoff of distance for lub

;; periodic systems
(define periodic   #t)      ; #f => non periodic, #t => periodic
(define ewald-tr   4.1)     ; time ratio Tr/Tk for Ewald summation
(define ewald-eps  1.0e-12) ; cut-off limit for Ewald summation
(define lattice    '(5.0  5.0  5.0)) ; size of the periodic box

;; iterative solver (ignored if flag-mat is true)
(define IT-solver "otmk") ; solver for libiter
; the following solvers are available
; "cg"       conjugate gradient method
; "cgs"      conjugate gradient squared (Weiss' Algorithm 11)
; "bicgstab" bi-conjugate gradient stabilized (Weiss' Algorithm 12)
; "sta"      bi-cgstab method
; "sta2"     bi-cgstab2 method
; "gpb"      gpbi-cg method
; "otmk"     orthomin method
; "gmres"    generalized minimum residual method
(define IT-max 2000)   ; max number of iterations
(define IT-n   20)     ; restart number
(define IT-eps 1.0e-6) ; accuracy of the solution
(define IT-debug 0)    ; set 1 to print the debug info for solve_iter()

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
(define a '())
; otherwise, poly codes are used in the calculations
;(define a #(
;1.0
;1.0
;1.0
;1.0
;1.0
;1.0
;1.0
;1.0
;))

; slip length
; set '() for no-slip system (slip = 0 for all particles)
(define slip '())
; otherwise, slip codes are used in the calculations
;(define slip #(
;1.0
;1.0
;1.0
;1.0
;1.0
;1.0
;1.0
;1.0
;))

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

; imposed-flow parameters
(define Ui '(0.0  0.0  0.0)) ; imposed translational velocity
(define Oi '(0.0  0.0  0.0)) ; imposed angular velocity relative to O
(define Ei '(0.0  0.0  0.0  0.0  0.0)) ; imposed strain relative to O

; auxiliary imposed-flow parameters for simple shear
(define shear-mode 0); 0 == imposed flow is given by Ui, Oi, Ei
                     ; 1 == simple shear (x = flow dir, y = grad dir)
                     ; 2 == simple shear (x = flow dir, z = grad dir)
(define shear-rate 0.0); the shear rate for shear-mode = 1 or 2
(define shear-shift 0.0); the initial condition of cell-shift
                        ; for shear-mode = 1 or 2

(define F0 '(0.0  0.0 -0.1)) ; applied force
(define T0 '(0.0  0.0  0.0)) ; applied torque

(define stokes 0.0)     ; effective stokes number
(define ncol   10)      ; frequency of collision check in dt for st != 0

;; Brownian dynamics parameters
(define peclet -1)      ; peclet number (negative means no Brownian force)
(define length 1.0)     ; unit of the length scale [micro m].
                        ; the parameters above are recognized by this.
                        ; just ignored for non-Brownian case (peclet<0).
(define BD-seed 0)      ; seed for random number generator for BD.

(define n-cheb-minv 50)      ; number of chebyshev coefficients for minv
(define n-cheb-lub  70)      ; number of chebyshev coefficients for lub

; time-integration scheme
(define BD-scheme "mid-point")
; the following algorithms are available
; "mid-point"        the mid-point algorithm
; "BanchioBrady03"   Banchio-Brady (2003)
; "BallMelrose97"    Ball-Melrose (1997)
; "JendrejackEtal00" Jendrejack et al (2000)
; "semi-implicit-PC" semi-implicit predictor-corrector
(define BB-n   100)    ; step parameter for Banchio-Brady-2003 algorithm

; dt-adjustment parameters
; NOTE: if 'rmin' above is defined by non-zero, dt-adjustment is just skipped.
(define BD-rmin 1.0)   ; overlap-param for dt-adjustment process in BD.
                       ; the condition is (r2 <= rmin * a2).
(define dt-lim 1.e-12) ; lower bound to shrink dt to prevent overlaps
                       ; set "dt" (or larger value) if you don't want 
                       ; to adjust "dt" but just reject it.

;; bond parameters
(define bonds '())
(define flag-relax #f) ; #f => stokesian dynamics with bond interactions
                       ; #t => relaxation dynamics only with bond interactions
(define gamma 1.0)     ; friction coefficient for relaxation dynamics

;; excluded volume parameters
(define ev-v   '())    ; v [(micro m)^3] for each chain type
(define ev-lim 5.0)    ; max distance for F^{EV} [micro m]
