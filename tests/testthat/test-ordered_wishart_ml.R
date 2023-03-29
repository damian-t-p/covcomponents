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

lme_mats <- get_covs(lme1)
