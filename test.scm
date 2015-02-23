
(use srfi-42)

(define-method covariance ((xs <vector>) (ys <vector>))  
  (let ([xa (mean xs)] [ya (mean ys)])
    (/ (sum (vector-map (^ (x y) (* (- x xa) (- y ya))) xs ys))
       (- (vector-length xs) 1))))

(define-method sum ((xs <vector>))
  (sum-ec (: i 0 (vector-length xs)) (~ xs i)))

(define-method length ((xs <vector>)) (vector-length xs))

(define-method mean ((xs <vector>)) (/ (sum xs) (length xs)))

(define-method + ((a <vector>) (b <vector>) . more)
  (apply + (vector-map + a b) more))

;; (let ([xs (list->vector (read-nums "xx2.txt"))])
(let ([xs (list->vector (iota 200000))])
  #?=(covariance xs xs))




