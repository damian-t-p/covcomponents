#' Computes REML or ML estimates of nested covariance components
#'
#' @param sos List of sum-of-squares matrices
#' @param constraints Matrix of order constraints for the factors
#' @param data An object inheriting `nesteddata`
#' @param ... Other arguments to pass to `multiordered_wishart_ml`
#'
#' @return A list of estimated covariance components
#'
#' @export
fit_balanced_nesteddata <- function(sos, constraints, data, ...) {
  
  full_covs <- multiordered_wishart_ml(sos, constraints, ...)

  cov_comps   <- list()
  
  prev_factor <- NULL
  divisor     <- 1

  # Split full covariance estimates into covariance components
  # Due to the constraints, each element of cov_comps will be non-negative definite
  for (factor in attr(data, "factors")) {

    if (is.null(prev_factor)) {
      cov_comps[[factor]] <- full_covs[[factor]]
    } else {
      cov_comps[[factor]] <- (full_covs[[factor]] - full_covs[[prev_factor]]) / divisor
    }
    
    cov_comps 

    divisor <- divisor * attr(data, "n_levels")[[factor]]
    prev_factor <- factor
    
  }

  cov_comps
}

#' Produces a list of identity matrices to use as initial guesses for covariance components
#'
#' @param data An object inheriting `nesteddata`
#'
#' @return A named list of identity matrices with an entry corresponding to each factor in `data`
init_covs_nesteddata <- function(data) {

  identity <- diag(x = 1, nrow = dim(data))

  init_covs <- rep(list(identity), attr(data, "n_factors"))
  names(init_covs) <- attr(data, "factors")

  init_covs
  
}
