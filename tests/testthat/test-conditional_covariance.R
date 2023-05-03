set.seed(1234)

df <- data.frame(
  dam  = paste0("d", c(1, 1, 1, 2, 2, 3, 3, 3, 3)),
  ind  = paste0("i", 1:9),
  t1   = rnorm(9),
  t2   = rnorm(9),
  t3   = rnorm(9)
)

data <- nesteddata(df, factors = c("ind", "dam"))

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
