

;; find argument min of y = f(x)
(define-method find-min-by-linear ((fun <procedure>) x step)
  (define (find-min x step min)
    (let ([v (fun (+ x step))])
      (if (> min v) (find-min 
  (define (next x direction)
    (let* ([opt-step (find-min x direction (fun x))]
           [X (+ x opt-step)])
      (if (>= step (abs opt-step)) X
          (next X (- direction)))))
  (next step step))

(define (calc-step-size next-gradient-next curr-gradient next curr)
  (let ([dg (- next-gradient curr-gradient)]
        [dx (- next curr)])
    (/ (inner-product dg dx)
       (inner-product dg dg))))

(find-min-by-linear (^ (x) (* x x)) 2 -0.1)

(define-method find-max-by-linear ((fun <procedure>) x)
  (find-min-by-linear (^ (t) (- (fun t))) x))