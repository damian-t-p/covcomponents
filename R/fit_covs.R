#' Estimate covariances for a covarianace components model
#'
#' @export
fit_covs <- function(data, ...) {
  UseMethod("fit_covs", data)
}

fit_covs.nesteddata <- function(data,
                                fixed_intercept = FALSE,
                                ...) {

  sos         <- sumofsquares(data, fixed_intercept)
  constraints <- linear_order(attr(data, "factors"))

  full_covs <- multiordered_wishart_ml(sos, constraints, ...)

  cov_comps   <- list()
  
  prev_factor <- NULL
  divisor     <- 1
  
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
