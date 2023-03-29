#' Generate pairwise constraints representing a linear order
#'
#' @param v A character vector of variable names in increasing order.
#'
#' @export
linear_order <- function(v) {

  n <- length(v)

  v_doubled = rep(v, each = 2)[-c(1, 2*n)]

  matrix(v_doubled, nrow = 2)
  
}
