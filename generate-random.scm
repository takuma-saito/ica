;; You can generate time series which obey AR(1) model.

(use math.matrix)
(use math.stat)
(use math.const)
(use srfi-27) ;; random
(use srfi-42) ;; loop

(define rand-mt
  (let1 r (make-random-source)
        (random-source-randomize! r)
        r))

;; 0.0 ~ 1.0
(define (rand-uniform)
  ((random-source-make-reals rand-mt)))

;; mean: 0, variance: 1.0
(define (rand-gauss-normal)
  (* (sqrt (* -2 (log (rand-uniform)))) (cos (* 2 pi (rand-uniform)))))

(define (rand-gauss mean variance)
  (+ mean (* variance (rand-gauss-normal))))

(define (rand-laplace u b)
  (let ([v (- 0.5 (rand-uniform))])
    (- u (* b (sgn v) (log (- 1.0 (* 2.0 (abs v))))))))

;; (rand-laplace3 1.0)

(define (make-series inits parameters variance)
  (define vars inits)
  (define (next)
    (let ((result ((apply$ +)
                   (map * (append '(1.0 1.0) vars)
                        (cons (rand-laplace 0.0 2.0) parameters)))))
      (set! vars (append (cdr vars) (list result)))
      result))
  next)

(define (generator series n)
  (list->vector
   (let loop ([counter n])
     (if (= counter 0) '()
         (cons
          (series) (loop (- counter 1)))))))

;; '#(#(2.4 3.3) #(-1.7 2.8))
(define in1 (generator (^ () (rand-laplace 0.0 1.0)) 1000))
(define in2 (generator (^ () (rand-laplace 0.0 1.0)) 1000))
(define m (vector->matrix '#(#(2.4 3.3) #(3.3 2.8))))
(define result (transpose (map (^ (a b) (* m (vector a b))) in1 in2)))
(define out1 (result 0))
(define out2 (result 1))
(define inv (inverse m))

;; debug
;; (normalize (m 0))
;; (normalize (m 1))
;; (p (normalize ((transpose inv) 1)) :print)
;; (normalize (inv 0))
;; (normalize (inv 1))
;; (p m :wolfram)
;; (p inv :print)

(define (out vars comment)
  (let ([xs (transpose vars)])
    (do-ec (: x xs)
           (print #`",(x 0) ,(x 1)"))
    (print "")))

(print #`"# input")
(print #`"# linear transform: ,(p m)")
(out (vector in1 in2) "# in")

(print #`"# output")
(print #`"# linear transform: ,(p m)")
(out (vector out1 out2) "# out")


