; test code for libstokes
; Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
; $Id: test-stokes.scm,v 1.1 2006/10/03 21:38:11 ichiki Exp $
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

(define sys (stokes-init))

(define np 8)
(define nm 8)
(stokes-set-np sys np nm)

(define lx 10.0)
(define ly 10.0)
(define lz 10.0)
(stokes-set-ll sys lx ly lz)

(define tratio 60.25)
(define zeta (zeta-by-tratio sys tratio))

(define cutlim 1.0e-12)
(stokes-set-zeta sys zeta cutlim)

(display "zeta = ")
(display zeta)
(newline)

;(stokes-lubcut-set sys 2.0000000001)
;(stokes-lubcut-get sys)
; with '-emit-setter' on swig, you can write those as
(set! (stokes-lubcut sys) 2.0000000001)
(stokes-lubcut sys)

(stokes-it-set sys (iter-init "gmres" 2000 20 1.0e-6 1))

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
(do ((i 0 (1+ i)))
    ((>= i np))
  (display i)
  (display " ")
  (display (darray-getitem pos (* i 3)))
  (display " ")
  (display (darray-getitem pos (+ 1 (* i 3))))
  (display " ")
  (display (darray-getitem pos (+ 2 (* i 3))))
  (newline))

(display "u:")
(newline)
(do ((i 0 (1+ i)))
    ((>= i np))
  (display i)
  (display " ")
  (display (darray-getitem u (* i 3)))
  (display " ")
  (display (darray-getitem u (+ 1 (* i 3))))
  (display " ")
  (display (darray-getitem u (+ 2 (* i 3))))
  (newline))

(set! (stokes-pos sys) pos)

(calc-res-ewald-3f sys u f)

(display "f:")
(newline)
(do ((i 0 (1+ i)))
    ((>= i np))
  (display i)
  (display " ")
  (display (darray-getitem f (* i 3)))
  (display " ")
  (display (darray-getitem f (+ 1 (* i 3))))
  (display " ")
  (display (darray-getitem f (+ 2 (* i 3))))
  (newline))