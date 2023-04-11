#' Use eigen to compute a square root of the non-negative definite symmetric matrix A
#'
#' The resulting U may not be symmetric, but satisfies `U'U == A`
#'
#' @param A Non-negative definite symmetric matrix
#' @param compute_inverse Logical indicating if inverse should also be computed. If this is `TRUE`, `A`
#' must be positive definite.
#'
#' @return If `compute_inverse` is `FALSE`, a matrix with dimensions matching `A`. If it is `TRUE`, a
#' list with entries names `sqrt` and `invsqrt`, containing a square root and inverse square root of
#' `A` respectively.
eigensqrt <- function(A, compute_inverse = FALSE) {
  decomp <- eigen(A, symmetric = TRUE)

  V <- decomp$vectors
  d <- decomp$values

  stopifnot(all(d > -1e-6))
  
  if (compute_inverse) {
    list(
      sqrt    = t(V) * sqrt(abs(d)),
      invsqrt = t(t(V) * 1 / sqrt(abs(d)))
    )
  } else {
    t(V) * sqrt(abs(d))
  }

}

#' Compute (A^(-1) + n E^(-1))^(-1) for a vector of ns
#'
#' @param A Symmetric nonnegative-definite square matrix
#' @param E Symmetric positive-definite square matrix
#' @param ns A vector of numerics
#' @param E_type If this equals `"prec"`, the `E` argument is instead `E^(-1)`
#'
#' @return A list of inverted matrices indexed by the vector `ns`
paired_inverse <- function(A, E, E_type = c("cov", "prec")) {

  E_type <- match.arg(E_type)

  U <- eigensqrt(A)

  if(E_type == "cov") {
    W <- U %*% solve(E, t(U))
  } else {
    W <- U %*% E %*% t(U)
  }

  W_eigen <- eigen(W, symmetric = TRUE)
  
  UtP <- t(U) %*% W_eigen$vectors

  inv_mats <- list()
  for(n in ns) {
    inv_mats[[as_label(n)]] <- UtP %*% diag(1/(W_eigen$values * n + 1)) %*% t(UtP)
  }

  inv_mats
  
}
