#' Compute sum-of-squares matrices from nested data
#'
#' Computes a list of between-factor sum-of-squares matrices (A[1], A[2], ..., A[n]).
#' If the model is correct, these are independent and have a Wishart distribution with
#'
#' A[k] ~ W(degf = (I[k] - 1) I[k] ... I[n],
#'          cov  = Sigma[1] + I[1] Sigma[2] + ... + I[1] ... I[n-1] Sigma[n]))
#'
#' @param data An object of class `nesteddata`
#' @param fixed_intercept Logical indicating whether the global mean in the model is a fixed effect.
#' This should be `TRUE` for maximum likelihood estimation. If `TRUE`, `A[n]` will have an additional
#' degree of freedom.
#'
#' @return A list of matrices indexed by factor name. Each matrix has a `degf` attribute, which encodes
#' the degrees of freedom in the corresponding likelihood function. This is usually the degrees of
#' freedom of the Wishart distribution, except in the case of ML estimation.
#'
#' @export
sumofsquares <- function(data,
                         fixed_intercept = FALSE) {

  stopifnot(is_balanced(data))
  
  noncentral_sos <- list()
  noncentral_sos[[attr(data, "factors")[1]]] <- data$sos_matrix

  central_mom <- list()
  
  ancestors <- rownames(data$group_sums)
  
  for (factor_idx in 2L:attr(data, "n_factors")) {
    factor       <- attr(data, "factors")[factor_idx]
    child_factor <- attr(data, "factors")[factor_idx - 1]

    children_per_level <- attr(data, "n_levels")[[child_factor]]
    
    # Compute the grouped means and corresponding sum-of-squares matrix
    means <- rowsum(data$group_sums, ancestors) / weight(attr(data, "n_levels"), factor_idx)
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

  # `factor` is now the highest level
  
  grand_mean <- colMeans(means)
  grand_sos  <- grand_mean %o% grand_mean
  
  children_per_level <- attr(data, "n_levels")[[factor]]
  
  central_mom[[factor]] <- structure(
    weight(attr(data, "n_levels"), factor_idx) *
      (noncentral_sos[[factor]] - grand_sos * children_per_level),
    degf = degf(attr(data, "n_levels"), factor_idx, fixed_intercept = fixed_intercept)
  )

  structure(
    central_mom,
    factors = attr(data, "factors")
  )
}

#' Degrees of freedom for nested designs
#'
#' Computes I[k] ... I[n-1] (I[n] - 1), or I[n] if k = n in the fixed-intercept case.
#'
#' @param ns Ordered vector of number of levels per factor
#' @param factor_idx Index of factor for which to compute degrees of freedom
#' @param fixed_intercept Logical indicating whether the model includes a fixed intercept
degf <- function(ns, factor_idx, fixed_intercept = FALSE) {
  ns <- ns[factor_idx:length(ns)]
  ns[1L] <- ns[1L] - !(fixed_intercept & (length(ns) == 1))
  prod(ns)
}

#' Number of individuals per level at a given factor
#'
#' Computes (I[1] ... I[k-1])
#'
#' @param ns Ordered vector of number of levels per factor
#' @param factor_idx Index of factor for which to compute weight
weight <- function(ns, factor_idx) {
  if (factor_idx == 1L) {
    1
  } else {
    prod(ns[1L:(factor_idx - 1L)])
  }
}
