#' ML estimates for Wishart covariance under linear order
#'
#' Given realisations of independent Wishart matrices A1 and A2, compute maximum likelihood estimates of
#' their underlying covariances Sigma1 and Sigma2 under the restriction that 0 <= Sigma1 <= Sigma2 with
#' respect to the non-negative definite ordering.
#'
#' This is done using an algorithm due to Amemiya (1985).
#'
#' @param A1,A2 Square matrices with a `degf` attribute.
#'
#' @return A list with entries `primal` and `dual`
#'
#' @export
ordered_wishart_ml <- function(A1, A2) {

  n1 <- attr(A1, "degf")
  n2 <- attr(A2, "degf")
  
  # sample covariances
  S1 <- A1 / n1
  S2 <- A2 / n2

  # S1 = U' U
  decomp <- eigen(S1, symmetric = TRUE)
  V <- decomp$vectors
  d <- decomp$values

  U     <- t(V) * sqrt(abs(d))
  U_inv <- t(t(V) * 1  /sqrt(abs(d)))

  # inv(U') S2 inv(U) = Q D Q'
  eigs <- eigen(t(U_inv) %*% S2 %*% U_inv, symmetric = TRUE)
  d <- eigs$values
  Q <- eigs$vectors

  # B S1 B' = I
  # B S2 B' = D
  B_inv <- t(U) %*% Q

  c_hat = 1 - n2/(n1 + n2) * pmax(1 - d, 0)
  d_hat = d + n1/(n1 + n2) * pmax(1 - d, 0)

  # primal solutions
  Sig1_hat <- B_inv %*% diag(c_hat, nrow = length(c_hat)) %*% t(B_inv)
  Sig2_hat <- B_inv %*% diag(d_hat, nrow = length(d_hat)) %*% t(B_inv)

  # dual solutions
  Psi1_hat <- (n1 / 2) * Sig1_hat - A1 / 2
  Psi2_hat <- (n2 / 2) * Sig2_hat - A2 / 2

  list(
    primal = list(Sig1_hat, Sig2_hat),
    dual   = list(Psi1_hat, Psi2_hat)
  )
  
}


#' ML estimates for Wishart covariance under general order
#'
#' Given a named list of realisations of independent Wishart matrices, compute the maximum likelihood
#' estimates of their underlying covariances under a given set of binary order constraints.
#'
#' This is done using an algorithm due to Calvin and Dykstra (1991)
#'
#' @param A_list A names list of square matrices with a `degf` attribute.
#' @param constraints A 2-row character matrix specifying the order constraints. Each column `v` encodes
#' a single order constraint through `v[1] <= v[2]`.
#' @param max.iter The number of EM iterations before termination.
#'
#' @export
multiordered_wishart_ml <- function(A_list,
                                    constraints,
                                    max.iter = 100) {

  p        <- nrow(A_list[[1]])
  zero_mat <- matrix(0, nrow = p, ncol = p)
  
  empty_dual      <- list(zero_mat, zero_mat)
  prev_iter_duals <- rep(list(empty_dual), ncol(constraints))

  # primal 
  prev_primal     <- list()
  primal          <- list()
  
  for (r in seq_len(max.iter)) {

    for (s in seq_len(ncol(constraints))) {

      constr <- constraints[, s]

      # Subtract previous dual for this constraint
      A_list[constr] <- mapply(
        \(A, Psi) {A - 2 * Psi},
        A_list[constr],
        prev_iter_duals[[s]],
        SIMPLIFY = FALSE
      )

      # Solve the current constraint only
      constr_sol <- ordered_wishart_ml(A_list[[constr[1]]], A_list[[constr[2]]])

      primal[constr]       <- constr_sol$primal
      prev_iter_duals[[s]] <- constr_sol$dual

      # Add current dual for this constraint
      A_list[constr] <- mapply(
        \(A, Psi) {A + 2 * Psi},
        A_list[constr],
        constr_sol$dual,
        SIMPLIFY = FALSE
      )
      
    }

    # check for convergence
    if (check_convergence(primal, prev_primal)) {
      break
    }

    prev_primal <- primal
    
  }

  primal
  
}

check_convergence <- function(mats1, mats2) {
  (min(length(mats1), length(mats2)) >= 1) &&
    all(mapply(\(m1, m2) all(abs(m1 - m2) < 1e-6), mats1, mats2))
}
