#' Uses the EM algorithm to maximise a given criterion in unbalanced data
#'
#' @param crit_argmax An function that maximises a criterion of interest given a list of arguments
#' @param arg_conditioner A function that computes the conditional expectation of sum-of-squares
#' matrices required for the `crit_argmax` function. Should take as its first argument the updating
#' #' parameter values.
#' @param init_params A list of initial parameter guesses
#' @param max.iter maximum iterations for the EM algorithm
#' @param check_convergence a function taking as arguments the current and previous parameters that
#' retugns a logical indicating whether the EM algorithm has converged.
#' @param crit_argmax_params List of non-updating arguments to pass to `crit_argmax`
#' @param sos_conditioner_params List of non-updating arguments to pass to `sos_conditioner`
#'
#' @return A list of parameters in the format of `init_params`
#' 
#' @export
unbalanced_EM <- function(crit_argmax,
                          sos_conditioner,
                          init_params,
                          max.iter               = 100L,
                          check_convergence      = \(curr_params, prev_params) FALSE,
                          crit_argmax_params     = list(),
                          sos_conditioner_params = list()) {

  prev_params <- init_params

  for (i in 1:max.iter) {

    cond_sos <- rlang::inject(
      sos_conditioner(prev_params, !!!sos_conditioner_params)
    )

    curr_params <- rlang::inject(
      crit_argmax(cond_sos, !!!crit_argmax_params)
    )

    if (check_convergence(curr_params, prev_params)) {
      break
    }

    prev_params <- curr_params
    
  }

  curr_params
  
}
