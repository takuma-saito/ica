(define-module math.matrix
  (use srfi-13) ;; string library
  (use srfi-42) ;; loop library
  (use math.std)
  (use srfi-43 :except (vector-map vector-map! vector-for-each vector-map-with-index))
  (export <matrix>
          eigens-sym transpose * - + / determinant
          solve decompose-LUP duplicate norm inverse
          vector->matrix make-init-vector
          power round-off othogonalize vec-plus
          map map$ map-with-index map-with-index$ normalize
          object-apply ref call call$ p vector->diag =: 
          for-each for-each-with-index inner-product eqv?))

(select-module math.matrix)

(define-class <matrix> ()
  ((m :init-keyword :m :init-value '#())))

(define-method ref ((A <matrix>) (i <integer>))
  (ref (ref A 'm) i))

(define-method ref ((A <matrix>) (i <integer>) (j <integer>))
  (ref (ref (ref A 'm) i) j))

(define-method object-apply ((A <matrix>) (i <integer>)) (ref A i))

(define-method object-apply ((A <matrix>) (i <integer>) (j <integer>)) (ref A i j))

(define-method (setter ref) ((obj <matrix>) (i <integer>) val)
  (vector-set! (ref obj 'm) i val))

(define-method (setter ref) ((obj <matrix>) (i <integer>) (j <integer>) val)
  (=: (ref A i j) val))

(define-method initialize ((A <matrix>) initargs)
  (next-method)
  (let-keywords
   initargs
   ([init 0] [row 0] [col 0] [m #f] . rest)
   (if (eqv? #f m)
       (=: (~ A 'm)
           ;; In gauche, you must initialize empty vector explicitly
           ;; to pass by value .
           (let ([vec (make-vector row)])
             (do-ec (: i 0 (vector-length vec))
                    (vector-set! vec i (make-vector col init)))
             vec)))))

;; (make <matrix> :col 5 :row 5 :init 0)

(define-method vector->matrix ((v <vector>)) (make <matrix> :m v))

(define-method map-with-index ((fun <procedure>) (A <matrix>))
    (vector-map-with-index
     (^ (i row)
        (vector-map-with-index
         (^ (j v)
            (fun i j v))
         row))
     (~ A 'm)))

(define-method map ((fun <procedure>) (A <matrix>))
  (make <matrix> :m (map-with-index (^ (i j v) (fun v)) A)))

(define-method eqv? ((A <matrix>) (B <matrix>))
  (every-with-index (^ (i j a) (eqv? (~ B i j) a)) A))

(define-method matrix->string
  ((A <matrix>) (col-sep <string>) (row-sep <string>))
  ($ (cut string-join <> row-sep) $
     vector->list
     (map (^ (x) ((.$ (cut string-join <> col-sep)
                      vector->list
                      (map$ number->string)) x))
          (~ A 'm))))

(define-method p ((A <matrix>) :optional (type :normal))
  (cond
   [(eqv? type :normal) (~ A 'm)]
   [(eqv? type :print)
    (print (matrix->string A "  " "\n"))]
   [(eqv? type :wolfram)
    ((cut (call$ string-concatenate) "{{" <> "}}")
     (matrix->string A ", " "}, {"))]))

(define-method for-each-with-index ((fun <procedure>) (A <matrix>))
  (begin (map-with-index fun A) '()))

(define-method for-each ((fun <procedure>) (A <matrix>))
  (begin (map fun A) '()))

(define-method for-each ((fun <procedure>) (v <vector>))
  (vector-for-each (^ (x) (fun x)) v))

(define-method every-with-index ((fun <procedure>) (A <matrix>))
  (call/cc
   (^ (return)
      (for-each-with-index
       (^ (i j x) (let ([out? (fun i j x)])
                (if (eqv? out? #f) (return #f))))
       A)
      #t)))

(define-method every ((fun <procedure>) (A <matrix>))
  (every-with-index (^ (i j x) (fun x)) A))

(define-method row-length ((A <matrix>)) (vector-length (~ A 'm)))

(define-method col-length ((A <matrix>)) (vector-length (~ A 0)))

(define-method row ((A <matrix>) (i <integer>)) (~ A i))

(define-method col ((A <matrix>) (i <integer>))
  (list->vector (list-ec (: j (row-length A)) (~ A j i))))

(define-method convolute-with-index  ((A <matrix>) (B <matrix>) fun)
  (make <matrix> :m (map-with-index (^ (i j v) (fun i j v (~ B i j))) A)))

(define-method convolute ((A <matrix>) (B <matrix>) fun)
  (make <matrix> :m (map-with-index (^ (i j v) (fun v (~ B i j))) A)))
  
(define-method + ((A <matrix>) (B <matrix>)) (convolute A B +))

(define-method - ((A <matrix>) (B <matrix>)) (convolute A B -))

(define-method * ((A <matrix>) (B <matrix>) . more)
  (apply *
         (convolute-with-index A B (^ (i j x y)
                                      (inner-product (row A i) (col B j))))
         more))

(define-method * ((A <matrix>) (n <number>)) (map (pa$ * n) A))

(define-method * ((n <number>) (A <matrix>)) (* A n))

(define-method * ((A <matrix>) (x <vector>))
  (list->vector
   (list-ec (: i 0 (vector-length x)) (inner-product (row A i) x))))

(define (min a b) (if (> a b) b a))

(define (make-unit-square-matrix n :optional (k n))
  (let [(m (make <matrix> :init 0 :row n :col k))]
    (do-ec (: i 0 (min k n)) (=: (~ m i i) 1))
    m))

(define-method swap! ((A <matrix>) (i <integer>) (j <integer>))
  (let ([a (~ A i)] [b (~ A j)])
    (vector-set! (~ A 'm) i b)
    (vector-set! (~ A 'm) j a)
    A))

(define-method duplicate ((A <matrix>))
  (let* ([m (~ A 'm)] [c (col-length A)] [r (row-length A)]
         [new-m (make <matrix> :init 0 :col c :row r)])    
    (do-ec (: i 0 r)
           (do-ec (: j 0 r)
                  (=: (~ new-m i j) (~ A i j))))
    new-m))

(define-method pivotize ((A <matrix>))
  (let ([M (make-unit-square-matrix (row-length A) (col-length A))])
    (do-ec (: i 0 (row-length A))
           (begin
             (let ([max -inf.0] [idx i])
               (do-ec (: j i (row-length A))
                      (if (> (A j i) max) (begin (=: max (A j i)) (=: idx j))))
               (if (not (eqv? idx i)) (swap! M idx i)))))
    M))

;; decompose A to PA = LU
(define-method decompose-LUP ((A <matrix>))
  (let* ([r (row-length A)] [c (col-length A)]
         [U (make <matrix> :init 0 :row r :col c)]
         [L (make <matrix> :init 0 :row r :col c)]
         [P (pivotize A)]
         [PA (* P A)])
    (do-ec (: i 0 (min r c)) (=: (~ L i i) 1))
    (do-ec (: i 0 r)
           (begin
             (do-ec (: j 0 i)
                    (=: (~ L i j)
                        ((cut / <> (U j j))
                         (- (~ PA i j)
                            ((apply$ +)
                             (list-ec (: k 0 j)
                                      (* (L i k) (U k j))))))))
             (do-ec (: j i c)
                    (=: (~ U i j)
                        (- (PA i j)
                           ((apply$ +)
                            (list-ec (: k 0 i)
                                     (* (L i k) (U k j)))))))))
    (values L U P)))

;; solve Ax = b, and return x
(define-method solve ((A <matrix>) (b <vector>))
  (let* ([r (row-length A)] [x (make-vector r 0)] [y (make-vector r 0)])
    (receive (L U P) (decompose-LUP A)
             (let ([new-b (* P b)])
               (do-ec (: i 0 r)
                      (=: (~ y i)
                          (- (new-b i)
                             ((apply$ +)
                              (list-ec (: k 0 i)
                                       (* (y k) (L i k)))))))
               (do-ec (: i (- r 1) -1 -1)
                      (=: (~ x i)
                          ((cut / <> (~ U i i))
                           (- (y i)
                              ((apply$ +)
                               (list-ec (: k i r)
                                        (* (x k) (U i k))))))))
               x))))

(define-method is-square? ((A <matrix>))
  (eqv? (row-length A) (col-length A)))

(define-method inverse ((A <matrix>))
  (if (not (is-square? A)) (error #`"square matrix is required")
      (let ([inv '()] [len (row-length A)])
        (do-ec (: i 0 len)
               (let ([vec (make-vector len 0)])
                 (=: (~ vec (- (- len 1) i)) 1)
                 (=: inv (cons (solve A vec) inv))))
        (make <matrix> :m (transpose (list->vector inv))))))

(define-method othogonalize ((x <vector>) basis-vectors)
  (- x
     ((apply$ +)
      (list-ec (: basis-vector basis-vectors)
               (* basis-vector
                  (/ (inner-product basis-vector x)
                     (inner-product basis-vector basis-vector)))))))

(define-method make-gramian-matrix ((x <vector>))
  (let* ([len (vector-length x)]
         [M (make <matrix> :init 0 :col len :row len)])
    (do-ec (: i 0 len)
           (do-ec (: j 0 len)
                  (=: (~ M i j) (* (~ x i) (~ x j)))))
    M))

;; (p (make-gramian-matrix '#(1 2 3 4 5)) :print)

(define-method determinant ((A <matrix>))
  (receive (L U P) (decompose-LUP A)
           (apply * (list-ec
                     (: i 0 (min (col-length A) (row-length A))) (U i i)))))

;; (* m7 (solve m7 '#(2.2 1.1 3.1)))
;; (* m3 (solve m3 '#(1 0 0)))
;; (* m4 (solve m4 '#(1 1 1)))

;; calculate ralyieght quotient
(define-method calc-rl-quotient ((A <matrix>) (x <vector>))
  (/ (inner-product (* A x) x) (inner-product x x)))

(define-method make-init-vector ((len <integer>))
  (let ([v (make-vector len 0)])
    (=: (~ v 0) 1)
    v))

;; power iteration method
;; return eigen-value, eigen-vector
(define-method minimum-eigen ((A <matrix>) :optional (init-value <integer>))
  (let* ([e +inf.0]
         [x (make-init-vector (row-length A))]
         [v init-value]
         [I (make-unit-square-matrix (row-length A))])
    (=: (~ x 0) 1)
    (while (> e 1.0e-12)
           (let ([new-x (solve (- A (* v I)) x)])
             (=: x (normalize new-x))
             (=: v (calc-rl-quotient A x))
             (=: e (norm (- (* A x) (* v x))))))
    (values v x)))

(define-method maximum-eigen ((A <matrix>))
  (let* ([e +inf.0]
         [x (make-init-vector (row-length A))]
         [v (calc-rl-quotient A x)])
    (while (> e 1.0e-12)
           (let ([new-x (* A x)])
             (=: x (normalize new-x))
             (=: v (calc-rl-quotient A x))
             (=: e (norm (- (* A x) (* v x))))))
    (values v x)))

(define-method transpose ((A <matrix>))
  (let* ([c (col-length A)]
         [r (row-length A)]
         [M (make <matrix> :init 0 :col (col-length A) :row (row-length A))])
    (do-ec (: i 0 r)
           (do-ec (: j 0 c)
                  (=: (~ M i j) (~ A j i))))
    M))

;; (transpose '((1 2) (3 4) (5 6)))

;; return eigen value and vector in ascent order.
(define-method eigens-n ((A <matrix>) (N <integer>))
  (let ([M (duplicate A)]
        [vecs '()]
        [vals '()])
    (do-ec (: i 0 N)
           (receive (val vec) (maximum-eigen M)
                    (=: M (- M (* val (make-gramian-matrix vec))))
                    (=: vals (cons val vals))
                    (=: vecs (cons vec vecs))))
    (values (list->vector vals) (list->vector vecs))))

(define-method is-sym? ((A <matrix>))
  (eqv? (transpose A) A))

;; return eingen values and eigen vector
(define-method eigens-sym ((A <matrix>))
  (if (not (is-sym? A)) (error "invalid matrix, you have to set symetric matrix")
      (receive (vals vecs) (eigens-n A (row-length A))
               (values vals (transpose (vector->matrix vecs))))))

(define-method vector->diag ((v <vector>))
  (let* ([len (vector-length v)]
         [M (make <matrix> :init 0 :col len :row len)])
    (do-ec (: i 0 len)
           (=: (~ M i i) (~ v i)))
    M))

;;
;; test
;;

(define v1 '#(1 2 3 4 5))
(define m0 (make <matrix> :init 0 :col 5 :row 5))
(define m1 (make <matrix> :m '#(#(1 2) #(3 4))))
(define m2 (make <matrix> :m '#(#(2 1) #(2 3))))
(define m3 (make <matrix> :m '#(#(12 -51 4) #(6 167 -68) #(-4 24 -41))))
(define m4 (make <matrix> :m '#(#(3 -9 8) #(8 2 -5) #(4 10 -9))))
(define m5 (make <matrix> :m '#(#(1 2 3) #(1 2 1) #(3 2 1))))
(define m6 (make <matrix> :m '#(#(0.1 0.8 0.9) #(2.3 0.3 5.5) #(2.2 3.4 3.5))))
(define m7 (make <matrix> :m '#(#(5 1 -2) #(1 6 -1) #(-2 -1 5))))
;; (p m7 :print)

;; (receive (vals Q) (eigens-sym m7)
;;          (p Q :print)
;;          (p (* Q (transpose Q)) :print)
;;          (p vals))
