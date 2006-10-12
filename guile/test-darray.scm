; test code for darray.scm
; Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
; $Id: test-darray.scm,v 1.1 2006/10/12 16:32:32 ichiki Exp $
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
(load "darray.scm")

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


;; set some darray 'pos'
(define np 8)
(define pos (new-darray (* np 3)))

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


(display "darray:")
(newline)
(display-darray3 pos np)


;; convert it to vector
(define vpos (darray->vector pos (* np 3)))

(newline)
(display "converted vector:")
(newline)
(display vpos)
(newline)

;; reconvert the vector to darray
(define pos2 (vector->darray vpos))

(newline)
(display "reconverted darray:")
(newline)
(display-darray3 pos2 np)
