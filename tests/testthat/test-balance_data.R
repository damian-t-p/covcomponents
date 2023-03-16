test_that("Name extension works", {
  expect_equal(extend_names(c(a = 2, b = 3, c = 0), existing_names = "a"),
               c(a.1 = "a", a.2 = "a", b = "b", b.1 = "b", b.2 = "b"))
})

set.seed(1234)

df <- data.frame(
  grp  = paste0("g", c(1, 1, 1, 1, 1, 2, 2, 2, 2)),
  sire = paste0("s", c(1, 1, 1, 2, 2, 3, 3, 3, 3)),
  dam  = paste0("d", c(1, 1, 2, 3, 4, 5, 6, 6, 7)),
  ind  = paste0("i", 1:9),
  t1   = rnorm(9),
  t2   = rnorm(9),
  t3   = rnorm(9)
)

data <- nesteddata(df, factors = c("ind", "dam", "sire", "grp"))
data_full <- add_unobs_levels(data)

means <- list(
  grp = matrix(
    c(100, 200), nrow = 2, ncol = 3,
    dimnames = list(c("g1","g2"))
  ),
  sire = matrix(
    c(10, 20, 30), nrow = 3, ncol = 3,
    dimnames = list(c("s1", "s2", "s3"))
  ),
  dam = matrix(
    c(1, 2, 3, 4, 5, 6, 7), nrow = 7, ncol = 3,
    dimnames = list(paste0("d", 1:7))
  )
)

balance(data, means)
