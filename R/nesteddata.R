#' Create a `nesteddata` object from a data frame
#'
#' @param df Wide data frame with one column for each level and one for each dimension of the
#' response variable.
#' @param factors Character vector containing the names of the level variables. Must be supplied in
#' the bottom-up order of nesting so that `factor[i]` is nested inside `factor[i+1]`.
#'
#' @return An object of class `nesteddata` containing the information necessary to compute covariance
#' component estimates.
#' The underlying object is a list with entries
#'   * `$sos_matrix` The matrix X'X, where X is the matrix of responses.
#'   * `$group_sums` Matrix of sums of responses grouped by the second-lowest factor.
#'
#' The object also contains the following attributes:
#'   * `attr(, "n_factors")` Integer recording the number of factors in the model. This will be 2 and
#' 3 for 1- and 2-way designs respectively.
#'   * `attr(, "factors")` Character vector of the names of the factors. The order determines the
#' nesting structure, with the factor `attr(, "factors")[i]` nested inside `attr(, "factors")[i+1]`.
#'   * `attr(, "parents")` List of character vectors indexed by the names of all factors except the
#' lowest and highest. The string `attr(, "parents")$factor["level"]` is the name of the level that
#' "level" is nested inside.
#'   * `attr(, "n_levels")` Integer vector recording the maximum number of different levels nested
#' in each parent for each factor.
#'   * `attr(, "n_observed")` List of integer vectors for every factor except the lowest. Records
#' the number of levels observed one level below each factor level.
#'
#' @export
nesteddata <- function(df, factors) {

  n_factors <- length(factors)
  stopifnot(n_factors >= 2L)
  
  value_matrix <- df %>%
    dplyr::select(-all_of(factors)) %>%
    as.matrix()

  label_matrix <- df %>%
    dplyr::transmute(across(all_of(factors), make.names)) %>%
    as.matrix()

  parents    <- list()
  n_observed <- list()
  
  for (i in 1L:(n_factors - 1L)) {
    curr_factor <- factors[i]
    
    curr_labels   <- label_matrix[, i]
    parent_labels <- label_matrix[, i + 1L]

    unique_curr_idx <- !duplicated(curr_labels)

    n_observed[[curr_factor]] <- c(tapply(curr_labels, parent_labels, \(v) length(unique(v))))
    
    if (i >= 2L) {
      parents[[curr_factor]] <- setNames(parent_labels[unique_curr_idx],
                                       curr_labels[unique_curr_idx])
    }
  }
  
  n_levels = c(
    sapply(n_observed, max),
    length(unique(label_matrix[, n_factors]))
  )
  names(n_levels) <- factors
  
  new_nesteddata(
    sos_matrix = t(value_matrix) %*% value_matrix,
    group_sums = rowsum(value_matrix, group = label_matrix[, 2L]),
    n_factors  = n_factors,
    factors    = factors,
    parents    = parents,
    n_levels   = n_levels, 
    n_observed = n_observed
  )
  
}

new_nesteddata <- function(sos_matrix,
                           group_sums,
                           n_factors,
                           factors,
                           parents,
                           n_levels,
                           n_observed) {

  stopifnot(is.integer(n_factors))  
  stopifnot(is.character(factors))
  stopifnot(is.list(parents) && all(sapply(parents, is.character)))
  stopifnot(is.integer(n_levels))
  stopifnot(is.list(n_observed) && all(sapply(n_observed, is.integer)))
  
  structure(
    list(sos_matrix = sos_matrix,
         group_sums = group_sums),
    class      = "nesteddata",
    n_factors   = n_factors,
    factors     = factors,
    parents    = parents,
    n_levels  = n_levels,
    n_observed = n_observed
  )
 
}
