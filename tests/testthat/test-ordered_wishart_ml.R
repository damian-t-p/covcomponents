set.seed(1234)


#one-way
p <- 3
df <- rnesteddata(
  covs = list(
    ind  = rnorm(p * p) %>% matrix(nrow = p) %>% {. %*% t(.)},
    dam  = rnorm(p * p) %>% matrix(nrow = p) %>% {. %*% t(.)}),
  nlevels_per_factor = c(ind = 2, dam = 3),
  var_names = c("t1", "t2", "t3")
)

data <- nesteddata(df, factors = c("ind", "dam"))
sos <- sumofsquares(data, method = "REML")

fit_covs.nesteddata(data, method = "REML")

ordered_wishart_ml(sos$ind, sos$dam)
halfsibdesign:::stepreml_1way(sos$ind, 3, sos$dam, 2)

set.seed(1234)
# two-way
p <- 3
df <- rnesteddata(
  covs = list(
    ind  = rnorm(p * p) %>% matrix(nrow = p) %>% {. %*% t(.)},
    dam  = rnorm(p * p) %>% matrix(nrow = p) %>% {. %*% t(.)},
    sire = rnorm(p * p) %>% matrix(nrow = p) %>% {. %*% t(.)}),
  nlevels_per_factor = c(ind = 2, dam = 2, sire = 3),
  var_names = c("t1", "t2", "t3")
)

data <- nesteddata(df, factors = c("ind", "dam", "sire"))

sos <- sumofsquares(data, method = "REML")

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
  
  resid_sds <- fit$sigma * sqrt(coef(fit$modelStruct$varStruct, FALSE, TRUE))
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

reml_mats <- fit_covs.nesteddata(data, method = "REML")
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


  1/2 * (
    n1 * (log(abs(det(solve(sig1) %*% G1))) - sum(diag(solve(sig1) %*% G1))) +
      n2 * (log(abs(det(solve(sig2) %*% G2))) - sum(diag(solve(sig2) %*% G2))) +
      n3 * (log(abs(det(solve(sig3) %*% G3))) - sum(diag(solve(sig3) %*% G3)))
  ) 

}

reml_score <- reml_crit(
  sos$ind / attr(sos$ind, "degf"),  reml_mats$ind,  2,
  sos$dam / attr(sos$dam, "degf"),  reml_mats$dam,  2,
  sos$sire / attr(sos$sire, "degf"), reml_mats$sire, 3
)

lme_score <- reml_crit(
  sos$ind / attr(sos$ind, "degf"),  lme_mats$ind,  2,
  sos$dam / attr(sos$dam, "degf"),  lme_mats$dam,  2,
  sos$sire / attr(sos$sire, "degf"), lme_mats$sire, 3
)

hsd_mats <- halfsibdesign::EM_fit(
  halfsibdesign::halfsibdata(df, df_format = "wide"),
  method = "REML")

hsd_score <- reml_crit(
  sos$ind / attr(sos$ind, "degf"),  hsd_mats$ind,  2,
  sos$dam / attr(sos$dam, "degf"),  hsd_mats$dam,  2,
  sos$sire / attr(sos$sire, "degf"), hsd_mats$sire, 3
)

test_that("Calvin-Dykstra produces non-negative definite estimates", {

  expect_true(all(eigen(reml_mats$ind)$values > -1e-6))
  expect_true(all(eigen(reml_mats$dam)$values > -1e-6))
  expect_true(all(eigen(reml_mats$sire)$values > -1e-6))
  
})

test_that("Calvin-Dykstra has a REML score no smaller than lme", {

  expect_true(reml_score >= lme_score - 1e-6)
  print(reml_score - lme_score)
  
})
