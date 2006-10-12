; convert utility between SCM Vector and darray in C
; Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
; $Id: darray.scm,v 1.1 2006/10/12 16:18:58 ichiki Exp $

; convert darray -> SCM Vector
(define (darray->vector x n)
  (define y (make-vector n))
  (do ((i 0 (1+ i)))
      ((>= i n))
    (vector-set! y i (darray-getitem x i)))
  y)

; convert SCM Vector -> darray
(define (vector->darray y)
  (define z (new-darray (vector-length y)))
  (do ((i 0 (1+ i)))
      ((>= i (vector-length y)))
    (darray-setitem z i (vector-ref y i)))
  z)
