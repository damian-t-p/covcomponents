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
                      group_effect = FALSE,
                      ...) {

  if (attr(data, "n_factors") == 2) {
    
    if (group_effect) {
      make_ccov_oneway_mainrandom(data, ...)
    } else {
      make_ccov_oneway_mainfixed(data, ...)
    }
  } else if (attr(data, "n_factors") == 3) {
    if(group_effect) {
      make_ccov_twoway_mainrandom(data, ...)
    } else {
      make_ccov_twoway_mainfixed(data, ...)
    }
    
  } else {
    stop("Conditional covariances only implemented for one- and two-way models",
         .call = FALSE)
  }
  
}

make_ccov_oneway_mainfixed <- function(data,
                                       prior_covs = identity_priors(data)) {

  factors <- attr(data, "factors")
  n_observed_unique <- unique(attr(data, "n_observed")[[ factors[1] ]])
  
  # $n[i]
  #  - the conditional variance of beta[i]
  ccov_blocks <- paired_inverse(
    prior_covs[[ factors[2] ]],
    prior_covs[[ factors[1] ]],
    n_observed_unique
  )

  n_lookup <- attr(data, "n_observed")$ind
  zero_mat <- zero_matrix(data)

  structure(
    function(rowname, colname) {
      if ((rowname %in% names(n_lookup)) && identical(rowname, colname)) {
        n_obs <- n_lookup[[rowname]]
        ccov_blocks[[as_label(n_obs)]]
      } else {
        zero_mat
      }
    },
    prior_covs = prior_covs
  )
  
}

make_ccov_oneway_mainrandom <- function(data,
                                        prior_covs = identity_priors(data),
                                        group_name = "group") {

  data_random_effect <- add_global_factor(data, group_name = group_name)

  make_ccov_twoway_mainfixed(
    data       = data_random_effect,
    prior_covs = c(prior_covs, rlang::list2({{group_name}} := flat_prior(data)))
  )
  
}

#' Create a list of unique numbers of observed levels per level two above a given index
promote_nobs <- function(data, level_idx = 1) {
  split(
    unname(sort(attr(data,"n_observed")[[ attr(data, "factors")[[level_idx]] ]])),
    attr(data,"parents")[[ attr(data, "factors")[[level_idx + 1]] ]]
  )
}

make_ccov_twoway_mainfixed <- function(data,
                                       prior_covs = identity_priors(data)) {

  factors <- attr(data, "factors")

  # number of observed individuals per middle factor
  n_observed_unique <- unique(attr(data, "n_observed")[[ factors[1] ]])

  prec_ind <- solve(prior_covs[[ factors[1] ]])

  unique_nobs_vecs <- unique(promote_nobs(data))
  
  # Compute conditional covariance of the alphas

  # List of the (Sigma_B + Sigma_E/n_j)^(-1) terms for all n_j
  W_summands <- paired_inverse(prec_ind,
                               prior_covs[[ factors[2] ]],
                               n_observed_unique,
                               E_type = "precision")

  W_summands <- mapply(`*`, W_summands, n_observed_unique, SIMPLIFY = FALSE)

  W <- list()
  for (nobs_vec in unique_nobs_vecs) {

    # Pick out the W_summands relevant to the current n_obs and summ
    W_sum <- Reduce(`+`, W_summands[paste0(nobs_vec)])

    W[[as_label(nobs_vec)]] <- paired_inverse(
      prior_covs[[ factors[3] ]],
      W_sum,
      1,
      E_type = "precision")[[1]]
    
  }
  
  
  # $n[i]
  #  - the conditional variance of beta[i]
  D_inv <- paired_inverse(prior_covs[[ factors[2] ]],
                          prior_covs[[ factors[1] ]],
                          n_observed_unique)

  W
  
}
