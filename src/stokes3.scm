; sample initialization file for stokes3
; SC lattice config of 8 particles in (5,5,5) box
; $Id: stokes3.scm,v 1.11 2008/04/26 19:06:25 kichiki Exp $

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
(define length 1.0)     ; the characteristic length [nm] (or [micro m]).
                        ; the parameters above are recognized by this.
                        ; the unit should be the same for parameters below
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
; an example of bonds
; note that the length unit should be the same for "length" above.
;(define bonds '(
;  (; bond 1
;   0         ; 1) spring type
;   (         ; 2) spring parameters (list with 3 elements)
;    0        ;    fene = 0 means (p1, p2) = (A^{sp}, L_{s})
;    1.0      ;    p1   = A^{sp}, scaled spring constant  (for fene == 0)
;    2.1)     ;    p2   = L_{s} / a, scaled max extension (for fene == 0)
;   ((0 1)    ; 3) list of pairs
;    (1 2)
;    (2 3))
;    -1)      ; 4) number of exclusion for lubrication
;             ;    negative means all particles in the chain is excluded.
;  (; bond 2
;   2         ; 1) spring type
;   (         ; 2) spring parameters (list with 3 elements)
;    1        ;    fene = 1 means (p1, p2) = (N_{K,s}, b_{K})
;    19.8     ;    p1 = N_{K,s}, the Kuhn steps for a spring (for fene = 1)
;    106.0)   ;    p2 = b_{K} [nm](or[micro m]), the Kuhn length (for fene = 1)
;   ((4 5)    ; 3) list of pairs
;    (5 6)
;    (6 7))
;     1)      ; 4) number of exclusion for lubrication
; ))
; where spring types are
;   0 : Hookean spring (Asp * (r - Ls)
;   1 : wormlike chain (WLC)
;   2 : inverse Langevin chain (ILC)
;   3 : Cohen's Pade approximation
;   4 : Warner spring
;   5 : Hookean spring (Asp * r / Ls)
;   6 : Hookean spring for dWLC

;; excluded volume parameters
; note that the length unit should be the same for "length" above.
(define ev-v   '())    ; v [nm^3] (or [micro m^3]) for each chain type
(define ev-lim 5.0)    ; max distance for F^{EV} [nm] (or [micro m])

;; angle parameters
(define angles '())
; an example
;(define angles '(
;  (; angle type 1
;   10.0    ; 1) constant (k^{angle})
;   0.0     ; 2) angle in degree (theta_0)
;   0       ; 3) scale flag (0 == scaled)
;           ;    in this case, the above value for k is just used.
;   ((0 1 2); 4) list of triplets
;    (1 2 3)
;    (2 3 4)
;   )
;  )
;  (; angle type 2
;   20.0    ; 1) constant (k^{angle})
;   90.0    ; 2) angle in degree (theta_0)
;   1       ; 3) scale flag (1 == not scaled yet)
;           ;    in this case, the potential is given by 
;           ;    (k/2) * kT * (theta - theta_0)^2
;   ((3 4 5); 4) list of triplets
;    (4 5 6)
;   )
;  )
;))

;; excluded volume in Debye-Huckel type
(define ev-dh   '())
; an example
; note that the length unit should be the same for "length" above.
;(define ev-dh '(
;  ; system parameters
;  4.0      ; 1) max distance for EV_DH interaction [nm] (or [micro m])
;  298.0    ; 2) temperature [K]
;  80.0     ; 3) dielectric constant of the solution
;  3.07     ; 4) Debye length [nm] (or [micro m])
;  (        ; 5) list of chain types
;   (; chain type 1
;    2.43    ; 1) nu [e/nm] (or [e/micro m])
;    5.00    ; 2) l0 [nm] (or [micro m])
;    (0 1 2) ; 3) list of particles
;   )
;   (; chain type 2
;    2.00    ; 1) nu [e/nm] (or [e/micro m])
;    4.00    ; 2) l0 [nm] (or [micro m])
;    (3 4)   ; 3) list of particles
;   )
;  )
;))
