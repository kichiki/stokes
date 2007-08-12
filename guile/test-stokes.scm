; test code for libstokes
; Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
; $Id: test-stokes.scm,v 1.8 2007/08/12 19:56:25 kichiki Exp $
;
; This program is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public License
; as published by the Free Software Foundation; either version 2
; of the License, or (at your option) any later version.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

(load-extension "./stokes.so" "SWIG_init")

; diaplay darray in 3
(define (display-darray3 da n)
  (do ((i 0 (1+ i)))
      ((>= i n))
    (display i)
    (display ": ")
    (display (darray-getitem da (* i 3)))
    (display " ")
    (display (darray-getitem da (+ 1 (* i 3))))
    (display " ")
    (display (darray-getitem da (+ 2 (* i 3))))
    (newline)))


(define sys (stokes-init))

(define np 8)
(define nm 8)
(stokes-set-np sys np nm)

(set! (stokes-periodic sys) 1) ; periodic boundary condition
(define lx 10.0)
(define ly 10.0)
(define lz 10.0)
(stokes-set-l sys lx ly lz)

(define ewald-tr 60.25)
(define xi (xi-by-tratio sys ewald-tr))

(define ewald-eps 1.0e-12)
(stokes-set-xi sys xi ewald-eps)

(display "xi = ")
(display xi)
(newline)

;(stokes-lubmin2-set sys 4.0000000001)
;(stokes-lubmin2-get sys)
; with '-emit-setter' on swig, you can write those as
(set! (stokes-lubmin2 sys) 4.0000000001)
;(stokes-lubmin2 sys)
(set! (stokes-lubmax sys) 4.0)
;(stokes-lubmax sys)

(stokes-set-iter sys "gmres" 2000 20 1.0e-6 1 (get-stdout))

(define pos (new-darray (* np 3)))
(define u   (new-darray (* np 3)))
(define f   (new-darray (* np 3)))

; set pos in SC lattice configuration
(darray-setitem pos 0  0.0) ; x component
(darray-setitem pos 1  0.0) ; y component
(darray-setitem pos 2  0.0) ; z component

(darray-setitem pos 3  5.0)
(darray-setitem pos 4  0.0)
(darray-setitem pos 5  0.0)

(darray-setitem pos 6  0.0)
(darray-setitem pos 7  5.0)
(darray-setitem pos 8  0.0)

(darray-setitem pos 9  0.0)
(darray-setitem pos 10 0.0)
(darray-setitem pos 11 5.0)

(darray-setitem pos 12 5.0)
(darray-setitem pos 13 5.0)
(darray-setitem pos 14 0.0)

(darray-setitem pos 15 0.0)
(darray-setitem pos 16 5.0)
(darray-setitem pos 17 5.0)

(darray-setitem pos 18 5.0)
(darray-setitem pos 19 0.0)
(darray-setitem pos 20 5.0)

(darray-setitem pos 21 5.0)
(darray-setitem pos 22 5.0)
(darray-setitem pos 23 5.0)

; set u and f
(do ((i 0 (1+ i)))
    ((>= i (* np 3)))
  (darray-setitem u i 1.0))

(display "pos:")
(newline)
(display-darray3 pos np)

(display "u:")
(newline)
(display-darray3 u np)

(stokes-set-pos sys pos)
(solve-res-3f sys u f)

;(define nc-f (stokes-nc-mob-f-init "test-stokes.res-3f.nc" np))
(define nc-f (stokes-nc-init "test-stokes.res-3f.nc"
			     np
			     0 ; nf
			     0 ; version
			     0 ; flag_poly
			     0 ; flag_Q
			     0 ; flag_it (time-dependent imposed flow)
			     ))
;; f0, x, u are active

(stokes-nc-set-f0 nc-f f)
(stokes-nc-set-time nc-f 0 0.0)
(stokes-nc-set-x nc-f 0 pos)
(stokes-nc-set-u nc-f 0 u)

(stokes-nc-free nc-f)


(display "f:")
(newline)
(display-darray3 f np)
