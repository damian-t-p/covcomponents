set.seed(1234)

df <- data.frame(
  sire = paste0("s", c(1, 1, 1, 2, 2, 3, 3, 3, 3)),
  dam  = paste0("d", c(1, 1, 2, 3, 4, 5, 6, 6, 7)),
  ind  = paste0("i", 1:9),
  t1   = rnorm(9),
  t2   = rnorm(9),
  t3   = rnorm(9)
)

data <- nesteddata(df, factors = c("ind", "dam", "sire"))

test_that("Factor structure inferred correctly", {
  expect_equal(attr(data, "n_factors"), 3) # 2-way model
  expect_mapequal(attr(data, "n_observed")$dam,
                  c(s1 = 2L, s2 = 2L, s3 = 3L))
  expect_mapequal(attr(data, "parents")$dam,
                  c(d1 = "s1", d2 = "s1", d3 = "s2", d4 = "s2", d5 = "s3", d6 = "s3", d7 = "s3"))
})
