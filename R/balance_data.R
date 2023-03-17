#' Use expected level means to balance an unbalanced design
#'
#' @param data An object of class `nesteddata`. This must only be unbalanced at the individual level,
#' with dummy levels present to represent unobserved levels.
#' @param level_means A list of matrices. The names of the list elements must be the same as those of
#' the group factors in data. Each list element should be a matrix whose row names are the levels of
#' the corresponding factor and whose rows are the relevant mean values.
#' @param global_mean A numeric vector. The global mean of the model.
#'
#' @return An object of class `nesteddata`, which is a balanced version of the input `data`.
#'
#' @export
balance <- function(data, level_means, global_mean = rep(0, dim(data))) {

  stopifnot(is_balanced(data, groups_only = TRUE))

  # Build conditional means for each dam
  # Since the unobserved individual component is always conditionally uncorrelated with the observed
  # data, the individual-level conditional means are always zero
  ancestors <- rownames(data$group_sums)
  ind_means <- 0 * data$group_sums + rep(global_mean, each = nrow(data$group_sums))
 
  for (factor in attr(data, "factors")[-1]) {
    ind_means <- ind_means + level_means[[factor]][ancestors, ]
    ancestors <- attr(data, "parents")[[factor]][ancestors]
  }

  # Compute conditional unobserved grouped sums by multiplying conditional means by number of unobserved
  # individuals
  lowest_factor <- attr(data, "factors")[1]
  n_missing <- attr(data, "n_levels")[lowest_factor] - attr(data, "n_observed")[[lowest_factor]]

  unobs_sums <- ind_means * n_missing

  # Now n_observed should indicate complete observations
  # Since data must already be group balanced, only individual-level balancing needs to be recorded
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
#' @param data An object of class `nesteddata`, which can have unbalanced groups factors.
#'
#' @return An object of class `nesteddata`, which has exactly the same information as the input `data`,
#' except with dummy levels with zero observations included where necessary to make the output
#' group-balanced.
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


#' Create a vector of new names that don't conflict with existing names
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
    names(new_names) <- make.names(full_names, unique = TRUE)[(k + 1):length(full_names)]
  }

  new_names
}
