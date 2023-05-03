cond_sos_oneway_mainfixed <- function(est_params, unbalanced_data) {

  ccov  <- make_ccov(unbalanced_data)
  cmean <- make_cmean(unbalanced_data,
                      ccov,
                      prior_global_mean = est_params$global_mean)

  balanced_data <- balance(data, cmean, cmean$global)

  sos <- sumofsquares(balanced_data,
                      method    = "ML",
                      intercept = TRUE)

  factors   <- attr(data, "factors")
  n_obs     <- attr(data, "n_observed")[[ factors[1] ]]
  n_max     <- attr(data, "n_levels")[[ factors[1] ]]
  top_max   <- attr(data, "n_levels")[[ factors[2] ]]
  n_missing <- n_max - n_obs
  
  top_levels <- names(n_obs)
  for (level in top_levels) {    
    sos[[ factors[1] ]] <- sos[[ factors[1] ]] +
      n_obs[[level]] * n_missing[[level]] / n_max *
      ccov(level, level)

    sos[[ factors[2] ]] <- sos[[ factors[2] ]] +
      n_missing[[level]]^2 * (1 / n_max) * (1 - 1 / top_max) *
      ccov(level, level)
      
  }
  
  sos[[ factors[1] ]] <- sos[[ factors[1] ]] +
    sum(n_missing) * (1 - 1 / n_max) *
    est_params[[ factors[1] ]]

  sos[[ factors[2] ]] <- sos[[ factors[2] ]] +
    sum(n_missing) * (1 / n_max) * (1 - 1 / top_max) *
    est_params[[ factors[1] ]]

  sos
  
}
