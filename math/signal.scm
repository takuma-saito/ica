(define-module math.signal
  (export read-series hanning-window make-series
          fcwt idft dft ifft fft reduce-size
          cut-bandwidth select-brandwidth))

(select-module math.signal)

(use math.const)
(use math.std)
(use srfi-42) ;; loop library
(use gauche.sequence)

(define (filter-with-index fn series)
  ((.$
    (map$ cadr)
    (filter$ (^ (x)
             (let ((index (car x)) (value (cadr x)))
               (fn index value)))))
   (zip (iota (length series) 0) series)))

(define (variance series)
  (let ((average (avg series)))
    (/ (sum (map (^ (x) (let ((v (- x average))) (* v v))) series))
       (length series))))

;; [max, min) :: range 
(define (make-series fn range step)
  (let* ((min (car range))
         (max (cadr range))
         (dt (/. (- max min) step)))
    (define (series n)
      (if (= n (- step 1)) (list (fn (* dt n)))
          (cons (fn (+ min (* dt n))) (series (+ n 1)))))
    (series 0)))

(define (filter-with-index fn series)
  ((map$ cadr)
   (filter (^ (x)
              (let ((index (car x)) (value (cadr x)))
                (fn index value)))
           (zip (iota (length series) 0) series))))

;; extract half items of array.
(define (reduce-size series)
  (let ((M (div (length series) 2)))
    (filter-with-index
     (^ (i x) (and (not (= i 0)) (> M i))) series)))

(define (pad-zero series)
  (let ((len (vector-length series)))
    (if (power2? len) series
        (let* ((num (- (ash (msb len) 1) len))
               (vector-zero (make-vector num)))
          (vector-fill! vector-zero 0)
          (vector-append series vector-zero)))))

;; 0.5*(1 + cos(2*pi*t/T))
(define (hanning x)
  (if (or (< x 0) (> x 1)) 0
      (* 0.5 (- 1 (cos (* 2.0 pi x))))))

(define (windowing window series)
  (let ((N (length series)))
    (map-with-index (^ (i x) (* x (window (/ i N)))) series)))

(define hanning-window
  (cut windowing hanning <>))

;; slow fourie transform
(define (dft series)
  (let* ([N (vector-length series)]
         [wave-fn (cut wave <> N -)])
    (map
     (^ (n)
        (sum
         (map-with-index
          (^ (k x) (* x (wave-fn (mod (* k n) N)))) series)))
     (vec-iota N))))

(define (idft series)
  (let* ([N (vector-length series)]
         [wave-fn (cut wave <> N +)])
     (map
      (^ (n)
         (sum
          (map-with-index
           (^ (k x) (/ (* x (wave-fn (mod (* k n) N))) N))
           series)))
      (vec-iota N))))

(define (merge even odd)
  (vector-map-with-index
   (^ (i _)
     (if (even? i)
         (vector-ref even (ash i -1))
         (vector-ref odd (ash i -1))))
   (make-vector (* 2 (vector-length even)))))

;; fast fourie transform
(define (fft-base fn series)
  (let* ((N (vector-length series))
         (h (ash N -1)))
    (if (= N 1) series
        (let ([even
               (fft-base fn
                         (vector-map-with-index
                          (^ (i _)
                            (+ (vector-ref series i)
                               (vector-ref series (+ i h))))
                          (make-vector h)))]
              [odd
               (fft-base fn
                         (vector-map-with-index
                          (^ (i _)
                            (* (fn i N)
                                (- (vector-ref series i)
                                   (vector-ref series (+ i h)))))
                          (make-vector h)))])
          (merge even odd)))))

(define (fft series)
  (fft-base (^ (n N) (wave n N -)) series))

(define (ifft series)
  (vector-map
   (cut / <> (vector-length series))
   (fft-base (^ (n N) (wave n N +)) series)))

;; auto-correlation function
(define (acf series)
  (let* ((N (vector-length series))
         (var (*. (variant series) N))
         (average (avg series))
         (s (normalize series)))
    (vector-map
     (^ (i _)
       ((cut / <> var)
        (sum
         (vector-map-with-index
          (^ (n x)
            (*. x (vector-ref s (mod (+ n i) N))))
          s))))
     s)))
    
(define (facf series)
  (let* ((s (normalize (pad-zero series)))
         (freq (fft s))
         (var (variant s))
         (N (vector-length s)))
    (vector-map
     (^ (x) (real-part (/ x var)))
     (ifft
      (vector->list
       (vector-product freq
                       (vector-map conjugate freq)))))))

;; todo: wavelet analysis

(define (morlet m z)
  (* (expt pi -0.25) (exp (* +i m z)) (exp (- (/ (* z z) 2.0)))))

(define (morlet-freq m z)
  (*. (expt pi -0.25)
      (if (> (real-part z) 0) 1 0)
      (expt e (/ (- (expt (- z m) 2)) 2))))

(define (cwt-wavelet wavelet series j k dt ds)
  (let* ((scale (* 2 dt (expt 2.0 (* j ds))))
         (coef (sqrt (/ dt scale))))
    ((cut / <> (vector-length series))
     (sum
      (vector-map-with-index
       (^ (i x)
         (* x
            coef
            (memo (^ (n j)
                    (conjugate
                     (wavelet (/ (* dt n) scale))))
                  (- i k) j)))
       series)))))

(define (with-scale fn series dt ds)
  (let* ((N (vector-length series))
         (K (x->integer (floor (/ (log (/ N 2)) (* ds (log 2.0)))))))
    (fn N K)))

(define (cwt-main wavelet printer series dt ds)
  (with-scale
   (^ (N K)
       (do-ec (: k 0 N)
             (do-ec (: j 0 K)
                   (printer (*. k dt) (*. j ds)
                            (cwt-wavelet wavelet series j k dt ds)))
             (print "")))
   series dt ds))

(define cwt
  (cut cwt-main <> <> <> <> 0.05))

(define (fcwt-wavelet wavelet-freq series scale dt)
  (let* ((freq (fft series))
         (N (vector-length series))
         (norm (* (sqrt (/ (* 2 pi scale) dt)) (/ 1 N))))
    (vector-map
     (^ (x) (* norm x))
     (ifft
      (vector-product
       (vector-map
        (^ (k)
          (let ((omega (/ (* 2 pi k) (* N dt))))
            (wavelet-freq (* scale omega))))
        (make-vector N))
       freq)))))

(define (fcwt-main wavelet-freq printer series dt ds)
  (with-scale
   (^ (N K)
     (let ((result (make-vector K)))
       (loop (k 0 K)
             (set! (ref result k)
                   (fcwt-wavelet wavelet-freq series
                                 (* 2 dt (expt 2.0 (* k ds))) dt)))
       (vector2d-print printer (vector-transpose result) dt ds)))
   series dt ds))

(define fcwt
  (cut fcwt-main <> <> <> <> 0.05))

;; (fcwt (cut morlet-freq 6.0 <>)
;;      (make-series (^ (x) (sin (* 2 5 pi x))) '(0 4) 16)
;;      (/. 4 16))

(define (pre-proc series)
  (hanning-window (normalize (pad-zero series))))

(define (phase-lag z)
  (atan (/ (real-part z) (imag-part z))))

;; generate normalized stream
(define (make-rand-norm)
  (define rand-real
    (random-source-make-reals
     (let1 r (make-random-source)
       (random-source-randomize! r)
       r)))
  (^ ()
    (* (sqrt (* -2.0 (log (rand-real))))
       (cos (* 2.0 pi (rand-real))))))

(define rand-norm (make-rand-norm))

(define (rand-norms n)
  (if (= n 0) '() (map (^ _ (rand-norm)) (iota n 0))))

(define (log10 x)
  (/ (log x) (log 10)))

(define (print-facf printer series)
  (let* ((result (facf (pre-proc series)))
         (rate (/. (vector-length result) (vector-length series))))
    (print-wave
     (^ (i x)
       (printer i x rate))
     (vector-map
      (^ (_ x) (real-part x)) result))))

(define (linear-fitting xs ys . opt)
  (let* ((norm-factor (if (null? opt) 0 (car opt)))
         (A (sum xs))
         (B (+ (+ (length xs) 1)
               (* norm-factor 2.0)))
         (C (sum ys))
         (D (+ (sum (map (^ (x) (* x x)) xs))
               (* norm-factor 2.0)))
         (E A)
         (F (sum (map * xs ys)))
         (a (/ (- (* C E) (* B F)) (- (* A E) (* D B))))
         (b (/ (- (* C D) (* A F)) (- (* B D) (* A E)))))
    (values a b)))

(define (cut-bandwidth L series)
  (let ((len (length series)))
    (map
     (^x (if (> (/ (* 2 (abs x)) len) L) x 0.0)) series)))

(define (select-brandwidth indexes series)
  (let ((len (- 1 (length series))))
    (map-with-index
     (^ (i x)
        (if (not (eqv? #f (find (^ (k)
                             (or (= i k)
                                 (= (- len i) k)))
                          indexes))) x 0.0))
     series)))

