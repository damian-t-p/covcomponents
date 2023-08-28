ml_oneway <- function(cond_sos, unbalanced_data) {

  sos <- multiordered_wishart_ml(cond_sos, linear_order(attr(unbalanced_data, "factors")))

  factors     <- attr(unbalanced_data, "factors")
  group_means <- unbalanced_data$group_sums / attr(unbalanced_data, "n_observed")[[ factors[1] ]]

  c(
    list(global_mean = colMeans(group_means)),
    sos
  )
  
}
