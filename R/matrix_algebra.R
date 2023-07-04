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

  nonzero <- as.numeric(d > 1e-6)
  diag_v <- sqrt(abs(d)) + (1 - nonzero)
  
  if (compute_inverse) {
    list(
      sqrt        = t(V) * sqrt(abs(diag_v)),
      invsqrt     = t(t(V) * 1 / sqrt(abs(diag_v))),
      nonzero_idx = nonzero
    )
  } else {
    t(V) * sqrt(abs(d))
  }

}

#' Compute (A^(-1) + n E^(-1))^(-1) for a vector of ns
#'
#' @param A Symmetric nonnegative-definite covariance matrix
#' @param E Symmetric positive-definite covariance or precision matrix
#' @param ns A vector of numerics
#'
#' @return A list of inverted matrices indexed by the vector `ns`
paired_inverse <- function(A, E, ns, keep_names = TRUE) {

  U <- eigensqrt(A)

  if (is_precision(E)) {
    W <- U %*% E %*% t(U)
  } else {
    W <- U %*% solve(E, t(U))
  }

  W_eigen <- eigen(W, symmetric = TRUE)
  
  UtP <- t(U) %*% W_eigen$vectors

  inv_mats <- list()
  for (n in ns) {
    inv_mat <- UtP %*% diag(1/(W_eigen$values * n + 1)) %*% t(UtP)

    if (keep_names) {
      dimnames(inv_mat) <- dimnames(A)
    }
    
    inv_mats[[as_label(n)]] <- inv_mat

  }

  inv_mats
  
}


covm <- function(mat) {
  structure(mat, type = "covariance")
}

precm <- function(mat) {
  structure(mat, type = "precision")
}

is_precision <- function(mat) {
  !is.null(attr(mat, "type")) && (attr(mat, "type") == "precision")
}
