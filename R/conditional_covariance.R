#' Compute the conditional covariances between effects in a nested model given observed data
#'
#' @param data An object of class `nesteddata`
#' @param group_effect A logical indicating wether the model includes an overall group random effect.
#' 
#' @return A closure `ccov(effect1, effect2)` taking specifications of effects. The effects should
#' be specified as character vectors with the nesting structure of the effect in question. For example,
#'
#' `ccov(c("sire1", "dam12"), "sire1")`
#'
#' returns the conditional covariance matrix `Cov(beta_12, alpha_1 | y^o)`. A `NULL` value of either
#' effect corresponds to the overall group effect.
make_ccov <- function(data,
                      group_effect = FALSE) {
  
}


make_ccov_list <- function(data) {
  
}
