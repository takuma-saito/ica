(define-module math.stat
  (use math.matrix)
  (use srfi-42)  ;; loop library
  (export sum length mean map covariance self-covariance self-cov cov
          average expected-value covariance-matrix cov-matrix sgn whitening
          corr correlation normalize-variable expected-vector
          whitening-matrix))

(select-module math.stat)

;; (sum-vec '#(#(1 2 3) #(2 2 2) (3 2 1))

