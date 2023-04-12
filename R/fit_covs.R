#' Estimate covariances for a covarianace components model
#'
#' @export
fit_covs <- function(data, ...) {
  UseMethod("fit_covs", data)
}

fit_covs.nesteddata <- function(data,
                                method    = c("ML", "REML"), 
                                intercept = TRUE,
                                ...) {

  method <- match.arg(method)

  if (method == "REML" & !intercept) {
    stop("The REML criterion assumes an intercept")
  }

  constraints <- linear_order(attr(data, "factors"))

  if (is_balanced(data)) {
    
    sos <- sumofsquares(data, method = method, intercept = intercept)
    fit_balanced_nesteddata(sos, constraints, data, ...)
    
  } else {

    unbalanced_EM(
      crit_argmax            = fit_balanced_nesteddata,
      crit_argmax_params     = list(constraints = constraints, data = data),
      arg_conditioner        = conditional_sos,
      arg_conditioner_parals = list(data = data),
      init_params            = init_covs_nesteddata(data),
      ...
    )
    
  }
  
  
}
