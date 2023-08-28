set.seed(1234)

df <- data.frame(
  dam  = paste0("d", c(1, 1, 1, 2, 2, 3, 3, 3, 3)),
  ind  = paste0("i", 1:9),
  t1   = rnorm(9),
  t2   = rnorm(9),
  t3   = rnorm(9)
)

data <- nesteddata(df, factors = c("ind", "dam"))

ests <- unbalanced_em(
  crit_argmax            = ml_oneway,
  sos_conditioner        = cond_sos_oneway_mainfixed,
  init_params            = c(identity_priors(data), list(global_mean = rep(0, 3))),
  crit_argmax_params     = list(unbalanced_data = data),
  sos_conditioner_params = list(unbalanced_data = data)
)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
