make_cmean <- function(data,
                       ccov,
                       prior_global_mean = zero_vector(data),
                       ...) {

  if (attr(data, "n_factors") == 2) {
    make_cmean_oneway_mainfixed(data, ccov, prior_global_mean, ...)
  }
  
}

make_cmean_oneway_mainfixed <- function(data,
                                        ccov,
                                        prior_global_mean = zero_vector(data)) {

  error_cov  <- attr(ccov, "prior_covs")[[attr(data, "factors")[1]]]
  error_prec <- precm(solve(error_cov))
  
  centered_sums <- data$group_sums - attr(data, "n_observed")[[attr(data, "factors")[1]]] %o% prior_global_mean
  skewed_sums   <- centered_sums %*% error_prec

  cond_mean <- skewed_sums * 0
  
  for (row in rownames(cond_mean)) {
    cond_mean[row, ] <- ccov(row, row) %*% skewed_sums[row, ]
  }

  top_name <- attr(data, "factors")[2]
  
  rlang::list2(
    global       = prior_global_mean,
    {{top_name}} := cond_mean
  )
  
}
