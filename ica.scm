;; This is independent analysis modules.

;; (use math.std)
(use math.matrix)
(use math.stat)
(use math.const)
(use srfi-42) ;; loop

(debug-print-width 500)

(define (radius->point x)
  (let ([z (expt e (* 0.0+i x))])
    (vector (real-part z) (imag-part z))))

(define (angle->point i)
  (radius->point (* 2 pi (/ i 360))))

(define (angle->radius i)
  (* (* 2 pi) (/ i 360)))

(define (point->radius point)
  (atan (/ (point 1) (point 0))))

(define (point->angle point)
  (* (/ (point->radius point) (* 2 pi)) 360))

(define (angle n fun)
  (do-ec (: i 0 n)
         (let ([v (fun (angle->point (* i (/ n 360))))])
           (print #`",(* i (/ n 360))   ,|v|"))))

(define angle-360 (cut angle 360 <>))
(define angle-3600 (cut angle 3600 <>))

;; kurt means kurtosis, zs must be whitenezed.
(define (kurt fun zs)
  (let ([trials (transpose zs)])
    (- (expected-value (^ (x) (let ([v (fun x)]) (* v v v v))) trials)
       3)))

(define (show-kurt converter vars)
  (let ([zs (whitening vars)])
    (converter
     (^ (weight)
        (kurt (^ (z) (inner-product z weight)) zs)))))

(define (show-gradient fun vars)
  (let ([zs (whitening vars)])
    (fun
     (^ (weight) (point->radius (normalize (calc-gradient weight zs)))))))

(define show-kurt-360 (cut show-kurt angle-360 <>))
(define show-kurt-3600 (cut show-kurt angle-3600 <>))
(define show-gradient-360 (cut show-gradient angle-360 <>))

(define (kurt-with-weight weight zs)
  (kurt (^ (z) (inner-product z weight)) zs))

(define (calc-gradient weight zs)
  (let ([trials (transpose zs)])
    (* 4 (- (expected-vector
             (^ (z) (* z (let ([v (inner-product z weight)]) (* v v v)))) trials)
            (* 3 weight)))))

(define (id x) x)

;; I use steepest decent method to find optimal weight.
(define (find-opt-weight zs :optional
                         (init-weight (make-init-vector (length zs)))
                         (otho id))
  (let* ([weight init-weight]
         [step-size 2.0]
         [curr-gradient (make-init-vector (length zs))])
    (do-ec (: i 0 10) ;; hard coded maximum loop
           (let* ([gradient (calc-gradient weight zs)]
                  [new-weight (otho (normalize (+ weight (* step-size gradient))))])
             ;; debug
             (display #`",(point->angle gradient) ,(normalize gradient)\n"
                      (standard-error-port))
             (flush-all-ports)
             (=: weight new-weight)))
    weight))

(define (transform V vars)
  (let ([trials (transpose vars)])
    (transpose
     (map (^ (var) (* V var)) trials))))

(define (shift items)
  (reverse (cdr (reverse items))))

(define (ica vars)
  (let ([weights '()]
        [zs (whitening vars)]
        [w (make-init-vector (length vars))])
    (do-ec (: i 0 (length vars))
           (begin
             (let ([new-weight
                    (find-opt-weight
                     zs w
                     (^ (weight) (othogonalize weight weights)))])
               (=: w (othogonalize new-weight (cons w weights)))
               (=: weights (cons new-weight weights)))))
    (vector->matrix
     (transpose (list->vector weights)))))

(define (is-comment? x)
  (eqv? (~ x 0) #\#))

(define (read-vars filename)
  (call-with-input-file
      filename
    (^ (port)
       (let* ([comments ""]
              [vars '()])
         (while (read-line port) (.$ not eof-object?) => line
                (cond [(string=? line "")] ;; continue
                      [(is-comment? line)
                       (=: comments (cons line comments))]
                      [else ;; Error sometimes occure in this line, future fix
                       (=: vars
                           (cons (list->vector
                                  (map string->number
                                       (string-split line " "))) vars))]))
         (values (transpose (list->vector (reverse vars)))
                 ((cut string-append <> "\n")
                  (string-join (reverse comments) "\n")))))))

(define (show-signals signals :optional (comments ""))
  (display comments)
  (do-ec (: signal (transpose signals))
         (print #`",(~ signal 0) ,(~ signal 1)")))

(define (with-estimate-signals filename fn)
  (receive (vars comments) (read-vars filename)
           (fn (transform (* (ica vars) (whitening-matrix vars)) vars)
               comments)))

(define (show-whitening filename)
  (receive (vars comments) (read-vars filename)
           (let ([zs (whitening vars)])
             (display comments)
             (do-ec (: trial (transpose zs))
                    (print #`",(trial 0)   ,(trial 1)")))))

(define (show-estimate-signals filename :optional (hook id))
  (with-estimate-signals filename
                         (^ (result comments)
                            (show-signals (map hook result) comments))))

(define (compare-estimate-signals outfile infile)
  (with-estimate-signals
   outfile
   (^ (signals comments)
    (do-ec (: sig signals)
           (let ([out '#()] [max -inf.0] [ans (read-vars infile)])
             (do-ec (: var ans)
                    (let ([c (abs (corr sig var))])
                      (if (> c max) (begin
                                      (=: out var)
                                      (=: max c)))))
             (let ([norm-sig (normalize-variable sig)]
                   [norm-out (normalize-variable out)])
               (show-signals (list norm-sig norm-out))
               (print "")))))))

;; (p (ica (read-vars "out.txt")) :print)

;; (show-estimate-signals "mixed.txt" normalize)

(compare-estimate-signals "out.txt" "in.txt")

;; debug
;; (show-whitening "mixed.txt")

;; (p (ica (read-vars "out.txt")) :print)
;; (point->angle '#(0.7116873719039802  0.718638469194167))
;; (point->angle '#(-0.7024963236006333  0.6953838872121385))
;; (show-signals (whitening (read-vars "out.txt")))

;; (show-kurt-360 (read-vars "out.txt"))
;; (show-gradient-360 (read-vars "out.txt"))
;; (show-whitening vars)


