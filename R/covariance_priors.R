zero_vector <- function(data) {
  setNames(rep(0, dim(data)), colnames(data$group_sums))
}

zero_matrix <- function(data) {
  zero_mat <- diag(0, dim(data))
  dimnames(zero_mat) <- dimnames(data$sos_matrix)

  zero_mat
}

identity_matrix <- function(data) {

  id_mat  <- diag(1, dim(data))
  dimnames(id_mat) <- dimnames(data$sos_matrix)
  
  id_mat
}

identity_priors <- function(data) {

  id_covs <- rep(list(identity_matrix(data)), attr(data, "n_factors"))
  names(id_covs) <- attr(data, "factors")

  id_covs

}

flat_prior <- function(data) {

  precm(zero_matrix(data))
  
}
