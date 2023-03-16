#' @export
balance <- function(data, level_means, global_mean = rep(0, dim(data))) {

  ancestors <- rownames(data$group_sums)
  ind_means <- 0 * data$group_sums + rep(global_mean, each = nrow(data$group_sums))
  
  for (factor in attr(data, "factors")[-1]) {
    ind_means <- ind_means + level_means[[factor]][ancestors, ]
    ancestors <- attr(data, "parents")[[factor]][ancestors]
  }

  lowest_factor <- attr(data, "factors")[1]
  n_missing <- attr(data, "n_levels")[lowest_factor] - attr(data, "n_observed")[[lowest_factor]]

  unobs_sums <- ind_means * n_missing

  balanced_n_obs <- attr(data, "n_observed")
  balanced_n_obs[[lowest_factor]][] <- attr(data, "n_levels")[lowest_factor]
  
  new_nesteddata(
    sos_matrix = data$sos_matrix + t(ind_means) %*% unobs_sums,
    group_sums = data$group_sums + unobs_sums,
    n_factors  = attr(data, "n_factors"),
    factors    = attr(data, "factors"),
    parents    = attr(data, "parents"),
    n_levels   = attr(data, "n_levels"),
    n_observed = balanced_n_obs
  )
}

#' Add empty unobserved factors to unbalanced data
#'
#' @export
add_unobs_levels <- function(data) {

  balanced_parents <- list()
  unbalanced_n_obs <- attr(data, "n_observed")
  
  # Iterate over all factors with children from highest to lowest
  for (factor_idx in (attr(data, "n_factors") - 1):2L) {

    factor       <- attr(data, "factors")[[factor_idx]]
    child_factor <- attr(data, "factors")[[factor_idx - 1]]
    
    new_names <- extend_names(
      attr(data, "n_levels")[factor] - unbalanced_n_obs[[factor]],
      names(attr(data, "parents")[[factor]])
    )

    balanced_parents[[factor]] <- c(attr(data, "parents")[[factor]], new_names)
    
    unbalanced_n_obs[[child_factor]] <- c(
      unbalanced_n_obs[[child_factor]],
      setNames(rep(0L, length(new_names)), names(new_names))
    )
  }

  # `factor` is now the lowest grouping factor

  # Add a row of empty sums for each newly-created lowest-grouping level
  group_sums <- rbind(
    data$group_sums,
    matrix(
      data     = 0,
      nrow     = length(new_names),
      ncol     = dim(data),
      dimnames = list(names(new_names))
    )
  )

  balanced_n_obs <- unbalanced_n_obs
  for (factor in attr(data, "factors")[-c(1L, attr(data, "n_factors"))]) {
    balanced_n_obs[[factor]][] <- attr(data, "n_levels")[[factor]]
  }

  new_nesteddata(
    sos_matrix = data$sos_matrix,
    group_sums = group_sums,
    n_factors  = attr(data, "n_factors"),
    factors    = attr(data, "factors"),
    parents    = balanced_parents,
    n_levels   = attr(data, "n_levels"),
    n_observed = balanced_n_obs
  )

}


#' Create new names
#'
#' @param n_missing Named integer vector detailing how many new names are to be created per parent.
#' @param existing_names Character vector of names that already exist, and so must be avoided.
extend_names <- function(n_missing, existing_names = character(0)) {
  new_names <- n_missing %>%
    mapply(\(count, parent) rep(parent, count), ., names(.)) %>%
    Reduce(c, .)

  k <- length(existing_names)
  full_names <- c(existing_names, new_names)

  if (length(new_names) >= 1L) {
    names(new_names) <- make.names(full_names, unique = TRUE)[(k+1):length(full_names)]
  }

  new_names
}
