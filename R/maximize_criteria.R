ml_oneway <- function(cond_sos, unbalanced_data) {

  sos <- multiordered_wishart_ml(cond_sos, linear_order(attr(unbalanced_data, "factors")))

  factors     <- attr(data, "factors")
  group_means <- data$group_sums / attr(data, "n_observed")[[ factors[1] ]]

  c(
    list(global_mean = colMeans(group_means)),
    sos
  )
  
}
