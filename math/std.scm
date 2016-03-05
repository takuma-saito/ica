(define-module math.std
  (use math.const)
  (use srfi-42)
  (export min max head =: map$ map-with-index$ sum
          sum-vec map-with-index transpose call call$
          object-apply eqv? inner-product + - * / norm
          normalize map p power round-off vec-iota id conjugate
          power2? msb wave))

(select-module math.std)

(define (min a b) (if (> a b) b a))

(define (max a b) (if (< a b) b a))

(define (wave i N type)
  (expt e (type (/ (*. 2.0 0.0+i i pi) N))))

(define head car)
(define tail cdr)
(define =: set!)

(define (power2? x)
  (= x (logand x (- x))))

(define (msb x)
  (if (power2? x) x
      (msb (logxor x (logand x (- x))))))

(define (conjugate z)
  (+ z (* 2.0 (- (real-part z) z))))

(define (id x) x)

(define (map$ fun) (^ (x) (map fun x)))

(define (map-with-index$ fun) (^ (x) (map-with-index fun x)))

(define =: set!)

(define-method sum ((xs <vector>))
  (sum-ec (: i 0 (vector-length xs)) (ref xs i)))

(define-method vec-iota ((h <integer>)) (list->vector (iota h 0)))

(define-method sum-vec ((xs <vector>))
  (let* ([len (vector-length xs)]
         [result (make-vector len 0)])
        (do-ec (: i 0 len)
               (=: result (vec-plus (ref xs i) result)))
        result))

(define-method map-with-index ((fun <procedure>) (A <vector>))
  (vector-map-with-index fun A))

(define-method transpose ((items <list>))
  (apply map (^ x x) items))

(define-method transpose ((items <vector>))
  (map list->vector
       (apply vector-map (^ x x) (vector->list items))))

(define-method object-apply ((A <vector>) (i <integer>)) (ref A i))

(define-method call ((fun <procedure>) . opts) (fun opts))

(define-method call$ ((fun <procedure>)) (^ opts (fun opts)))

(define-method eqv? ((A <vector>) (proc <procedure>))
  (every (^ (x) (proc x)) A))

(define-method inner-product ((v1 <vector>) (v2 <vector>))
  (sum-ec (: i 0 (vector-length v1)) (* (ref v1 i) (ref v2 i))))

(define-method + ((x <vector>) (y <vector>) . rest)
  (apply + (vector-map + x y) rest))

(define-method + ((a <vector>) (b <vector>) . more)
  (apply + (vector-map + a b) more))
(define-method - ((a <vector>) (b <vector>)) (vector-map - a b))
(define-method * ((a <vector>) (b <vector>) . more)
  (apply * (vector-map * a b) more))
(define-method / ((a <vector>) (b <number>))
  (vector-map (^ (v) (/ v b)) a))
(define-method - ((a <vector>) (b <number>))
  (vector-map (^ (v) (- v b)) a))
(define-method * ((items <vector>) (n <number>))
  (vector-map (^ (v) (* v n)) items))
(define-method * ((n <number>) (items <vector>))
  (vector-map (^ (v) (* v n)) items))
(define-method norm ((items <vector>))
  (sqrt (inner-product items items)))

(define-method normalize ((A <vector>))
  (/ A (norm A)))

(define-method length ((xs <vector>)) (vector-length xs))
(define-method mean ((xs <vector>)) (/ (sum xs) (length xs)))
(define average mean)
(define-method map ((fun <procedure>) (xs <vector>) (ys <vector>))
  (vector-map (^ (x y) (fun x y)) xs ys))

(define-method covariance ((xs <vector>) (ys <vector>))  
  (let ([xa (mean xs)] [ya (mean ys)])
    (/ (sum (map (^ (x y) (* (- x xa) (- y ya))) xs ys)) (- (length xs) 1))))

(define-method self-covariance ((xs <vector>))
  (covariance xs xs))

(define cov covariance)
(define self-cov self-covariance)

(define-method correlation ((xs <vector>) (ys <vector>))
  (/ (cov xs ys) (* (cov xs xs) (cov ys ys))))

(define corr correlation)

(define-method map ((fun <procedure>) (v <vector>)) (vector-map fun v))
(define-method map ((fun <generic>) (v <vector>)) (vector-map fun v))

(define-method p ((items <pair>) :optional (type :print))
  (define (to-s items)
    (string-join (map x->string items) "\n"))
  (cond [(eqv? type :string) (to-s items)]
        [(eqv? type :print) (display (to-s items))]))

(define-method p ((v <vector>) :optional opts)
  (p (vector->list v) opts))

(define-method power ((x <number>) (n <integer>))
  (let ([ret 1] [v x])
    (while (not (eqv? n 0))
           (if (eqv? (logand n 1) 1) (=: ret (* v ret)))
           (=: v (* v v))
           (=: n (ash n -1)))
    ret))

(define-method round-off ((x <number>) (n <integer>))
  (let ([N (power 10 n)])
    (/ (floor (* (inexact x) N)) N)))

(define-method print ((vec <vector>)))

(define-method linear-normalize ((vec <vector>))
  (let ([m (mean vec)]
        [var (sqrt (variance vec))])
    (map (.$ (cut / <> var) (cut - <> m)) vec)))

(define (expected-value fun var)
  (/ (sum (map fun var)) (length var)))
(define (expected-vector fun var)
  (/ (sum-vec
      (vector-map fun var))
     (vector-length var)))

;; vars is <vector> list, they are random variables
(define (covariance-matrix vars)
  (let* ([len (length vars)]
         [M (make <matrix> :init 0 :col len :row len)])
    (do-ec (: i 0 len)
           (do-ec (: j 0 len)
                    (=: (~ M i j) (cov (vars i) (vars j)))))
    M))

(define cov-matrix covariance-matrix)

(define-method sgn ((x <number>))
  (cond [(> x 0) 1]
        [(< x 0) -1]
        [(= x 0) 0]))

(define printer (cut p <> :print))

(define (whitening-matrix vars)
  (receive (vals Q) (eigens-sym (cov-matrix vars))
           ;; #?=vals (p Q :print)
           (* (vector->diag (map (^ (x) (/ 1.0 (sqrt x))) vals))
              (transpose Q))))

;; nomalizing random vars, mean = 1 and covariance-matrix = I
(define (whitening vars)
  (let ([V (whitening-matrix vars)] [mn (map mean vars)])
    ((.$ transpose
         (map$ (pa$ * V))
         (map$ (pa$ - mn))
         transpose)
     vars)))

;; normalize 1 variable mean = 0, sigma = 1
(define-method normalize-variable ((var <vector>))
  (let ([mn (mean var)]
        [sigma (sqrt (self-cov var))])
    (map (^ (x) (* (/ (- x mn) sigma))) var)))


;; separator is space
;; <string> => <vector>
(define-method read-nums ((filename <string>))
  (call-with-input-file
      filename
    (^ (port)
       (let ([ret '()])
         (while (read-line port) (.$ not eof-object?) => line
                (let ([xs ($ (map$ x->number) $
                             (filter$ (pa$ (.$ not equal?) ""))
                             (string-split line " "))])
                  (=: ret
                      (case (length xs)
                        [(0) ret] ;; skip
                        [(1) (cons (car xs) ret)]
                        [else (cons (list->vector xs) ret)]))))
         (list->vector (reverse ret))))))
  
