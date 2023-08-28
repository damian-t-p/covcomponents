#' Sample a random balanced nested design
#'
#' @param covs List of covariance components
#' @param nlevels_per_factor Named numeric vector indicating the number of distinct levels at each
#' level of nestine
#' @param gloal_mean Numeric vector. The fixed effect of the model
#' @param factors The names of the factors in increasing order of the nesting structure
#' @param level_names A list of level names for each factor
#' @param var_names A character vector of the names of the variables
#' @param format A character vector indicating the type of outpu
#'
#' @return Either a tibble or an object of class `nesteddata`, depending on the value of `format`.
#' 
#' @export
rnesteddata <- function(covs,
                        nlevels_per_factor,
                        global_mean = rep(0, nrow(covs[[1]])),
                        factors     = names(covs),
                        level_names = NULL,
                        var_names   = rownames(covs[[1]]),
                        format      = c("table", "data")) {

  format <- match.arg(format)

  if ( is.null(names(covs)) ) {
    names(covs) <- factors

    if ( is.null(factors) ) {
      stop(
        "covariance components must be names",
        call. = FALSE
      )
    }
    
  } else {
    
    if ( !identical(names(covs), names(nlevels_per_factor)) ) {
      stop(
        "covs and nlevels_per_factor must have the name structure",
        call. = FALSE
      )
    }

    if ( !setequal(factors, names(covs)) ) {
      stop(
        "factor format must have the same entries as the names of covs and nreps"
      )
    }
    
  }

  # Order parameters in increasing order of nesting
  covs <- covs[factors]
  nlevels_per_factor <- nlevels_per_factor[factors]

  nreps  <- rev(cumprod(rev(nlevels_per_factor)))
  total_reps <- nreps[[1L]]
  ndims <- nrow(covs[[1L]])
  
  # List of matrices comprising the levels of the data
  level_vals <- list()
  
  for (factor in factors) {
    white_data <- matrix(
      rnorm(nreps[[factor]] * ndims),
      nrow = nreps[[factor]],
      ncol = ndims
    )

    cov_data <- white_data %*% eigensqrt(covs[[factor]])

    if (is.null(level_names[[factor]])) {
      curr_names <- paste0(factor, seq_len(nreps[[factor]]))
    } else {
      curr_names <- level_names[[factor]]
    }

    rownames(cov_data) <- curr_names

    # Repeat each row so that the matrix has nreps rows
    expanded_names <- rep(curr_names, each = total_reps / nreps[[factor]])

    level_vals[[factor]] <- cov_data[expanded_names, ]
  }

  # Add all of the component matrices together
  full_vals <- Reduce(
    f    = `+`,
    x    = level_vals,
    init = matrix(global_mean, nrow = total_reps, ncol = ndims, byrow = TRUE)
  )

  if (is.null(var_names)) {
    var_names <- paste0("var", seq_len(ndims))
  }

  colnames(full_vals) <- var_names

  name_table <- tibble::as_tibble(lapply(level_vals, rownames))
  full_table <- cbind(name_table, tibble::as_tibble(full_vals))
  
  if (format == "table") {
    full_table
  } else {
    nesteddata(full_table, factors)
  }
  
}
