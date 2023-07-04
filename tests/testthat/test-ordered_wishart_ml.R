set.seed(1234)

df <- data.frame(
  sire = paste0("s", c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)),
  dam  = paste0("d", c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)),
  ind  = paste0("i", 1:12),
  t1   = rnorm(12),
  t2   = rnorm(12),
  t3   = rnorm(12)
)

data <- nesteddata(df, factors = c("ind", "dam", "sire"))
sos  <- sumofsquares(data)

ordered_wishart_ml(sos$ind, sos$dam)

cd_mats <- multiordered_wishart_ml(
  A_list      = sos,
  constraints = linear_order(c("ind", "dam", "sire"))
)

suppressWarnings({

  df_long <- tidyr::pivot_longer(
    df,
    dplyr::starts_with("t"),
    names_to = "trait",
    values_to = "value"
  )
  
  lme1 <- nlme::lme(
    fixed       = value ~ -1 + trait,
    data        = df_long,
    random      = ~ -1 + trait | sire/dam,
    weights     = nlme::varIdent(form = ~ 1 | trait),
    correlation = nlme::corSymm(form =  ~ 1 | sire/dam/ind),
    method      = "REML",
    control     = nlme::lmeControl(opt = "optim", returnObject = TRUE)
  )

})

get_covs <- function(fit){
  
  resid_sds <- fit$sigma * sqrt((1 + c(0, fit$modelStruct$varStruct)))
  resid_corr <- as.matrix(fit$modelStruct$corStruct)[[1]]
  S1 <- outer(resid_sds, resid_sds) * resid_corr
  
  S2 <- fit$modelStruct$reStruct$dam %>%
      as.matrix()  %>%
      {. * fit$sigma^2}
  rownames(S2) <- NULL
  colnames(S2) <- NULL
  
  S3 <- fit$modelStruct$reStruct$sire %>%
      as.matrix()  %>%
      {. * fit$sigma^2}
  rownames(S3) <- NULL
  colnames(S3) <- NULL
  
  out <- list(
    ind  = S1,
    dam  = S2,
    sire = S3
  )
  
  return(out)
}

reml_mats <- fit_covs(data, method = "REML")
lme_mats  <- get_covs(lme1)


reml_crit <- function(M1, S1, I1, M2, S2, I2, M3, S3, I3) {
  p <- nrow(M1)

  n1 <- I3*I2*(I1 - 1)
  n2 <- I3*(I2 - 1)
  n3 <- I3 - 1
  n <- n1 + n2 + n3

  G1 <- M1/n1
  G2 <- M2/n2
  G3 <- M3/n3

  sig1 <- S1
  sig2 <- S1 + I1 * S2
  sig3 <- S1 + I1 * S2 + I1*I2 * S3

  c <- -p*n/2 * log(2*pi) +
    -p/2 * log(I1*I2*I3) +
    #-1/2 * log(det(2*pi*sig3)) +
    -1/2 * (
      n1 * log(det(G1)) +
        n2 * log(det(G2)) +
        n3 * log(det(G3))
    )

  1/2 * (
    n1 * (log(abs(det(solve(sig1) %*% G1))) - sum(diag(solve(sig1) %*% G1))) +
      n2 * (log(abs(det(solve(sig2) %*% G2))) - sum(diag(solve(sig2) %*% G2))) +
      n3 * (log(abs(det(solve(sig3) %*% G3))) - sum(diag(solve(sig3) %*% G3)))
  ) + c

}

reml_score <- reml_crit(
  sos$ind / 6,  reml_mats$ind,  2,
  sos$dam / 3,  reml_mats$dam,  2,
  sos$sire / 2, reml_mats$sire, 3
)

lme_score <- reml_crit(
  sos$ind / 6,  lme_mats$ind,  2,
  sos$dam / 3,  lme_mats$dam,  2,
  sos$sire / 2, lme_mats$sire, 3
)

test_that("Calvin-Dykstra produces non-negative definite estimates", {

  expect_true(all(eigen(reml_mats$ind)$values > -1e-6))
  expect_true(all(eigen(reml_mats$dam)$values > -1e-6))
  expect_true(all(eigen(reml_mats$sire)$values > -1e-6))
  
})

test_that("Calvin-Dykstra has a REML score no smaller than lme", {

  expect_true(reml_score >= lme_score - 1e-2 * abs(lme_score))
  
})
