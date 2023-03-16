#' Compute sum-of-squares matrices from nested data
#'
#' Wishart A parameterisation
sumofsquares <- function(data,
                         fixed_intercept = FALSE) {

  noncentral_sos <- list()
  noncentral_sos[[attr(data, "factors")[1]]] <- data$sos_matrix

  central_mom <- list()
  
  ancestors <- rownames(data$group_sums)
  
  for (factor_idx in 2L:attr(data, "n_factors")) {
    factor <- attr(data, "factors")[factor_idx]
    child_factor <- attr(data, "factors")[factor_idx - 1]

    children_per_level <- attr(data, "n_levels")[[child_factor]]
    
    # Compute the grouped means and corresponding sum-of-squares matrix
    means <- rowsum(data$group_sums, ancestors) * weight(attr(data, "n_levels"), factor_idx)
    noncentral_sos[[factor]] <- t(means) %*% means

    # Centre the previous sum-of-squares matrix
    # cf (ind_sos - K dam_sos)
    central_mom[[child_factor]] <- structure(
      weight(attr(data, "n_levels"), factor_idx - 1) *
        (noncentral_sos[[child_factor]] - noncentral_sos[[factor]] * children_per_level),
      degf = degf(attr(data, "n_levels"), factor_idx - 1)
    )
    
    ancestors <- attr(data, "parents")[[factor]][ancestors]
  }

  grand_mean <- colMeans(means)
  grand_sos  <- grand_mean %o% grand_mean
  children_per_level <- attr(data, "n_levels")[[factor]]
  
  central_mom[[factor]] <- structure(
    weight(attr(data, "n_levels"), factor_idx) *
      (noncentral_sos[[child_factor]] - grand_sos * children_per_level),
    degf = degf(attr(data, "n_levels"), factor_idx, fixed_intercept = fixed_intercept)
  )

  structure(
    central_mom,
    factors = attr(data, "factors")
  )
}

#' Compute degrees of freedom
#'
#' @param ns Ordered vector of number of levels per factor
#' @param factor_idx Index of factor for which to compute degrees of freedom
#' @param fixed_intercept Logical indicating whether the model includes a fixed intercept
degf <- function(ns, factor_idx, fixed_intercept = FALSE) {
  ns <- ns[factor_idx:length(ns)]
  ns[1L] <- ns[1L] - !fixed_intercept
  prod(ns)
}

weight <- function(ns, factor_idx) {
  if (factor_idx == 1L) {
    1
  } else {
    1/prod(ns[1L:(factor_idx - 1L)])
  }
}
