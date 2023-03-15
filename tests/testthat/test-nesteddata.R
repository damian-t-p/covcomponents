set.seed(1234)

df <- data.frame(
  sire = c(1, 1, 1, 2, 2, 3, 3, 3, 3),
  dam  = c(1, 1, 2, 3, 4, 5, 6, 6, 7),
  ind  = 1:9,
  v1   = rnorm(9),
  v2   = rnorm(9),
  v3   = rnorm(9)
)

data <- nesteddata(df, levels = c("ind", "dam", "sire"))

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
