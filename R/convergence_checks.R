#' Generate a closure that checks for convergence of matrices in the l-infinity norm
#'
#' @param eps The threshold at which convergence is reaced
#' 
#' @export
converged_l_infinity <- function(eps = 1e-6) {
  function(mats1, mats2) {
    (min(length(mats1), length(mats2)) >= 1) &&
      all(mapply(\(m1, m2) all(abs(m1 - m2) < eps), mats1, mats2))
  }
}
