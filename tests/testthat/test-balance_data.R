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

test_that("Group balancing works", {
  expect_false(is_balanced(data))
  expect_false(is_balanced(data, groups_only = TRUE))
  expect_false(is_balanced(data_full))
  expect_true(is_balanced(data_full, groups_only = TRUE))
})

means <- list(
  grp = matrix(
    c(100, 200), nrow = 2, ncol = 3,
    dimnames = list(c("g1", "g2"))
  ),
  sire = matrix(
    c(10, 20, 30, 0), nrow = 4, ncol = 3,
    dimnames = list(c("s1", "s2", "s3", "g2"))
  ),
  dam = matrix(
    c(1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0), nrow = 12, ncol = 3,
    dimnames = list(c(paste0("d", 1:7), "s1", "s2", "g2", "g2.1", "g2.2"))
  )
)

data_balanced <- balance(data_full, means)

test_that("Individual balancing works", {
  expect_true(is_balanced(data_balanced))
  expect_true(is_balanced(data_balanced, groups_only = TRUE))
})

test_that("Balancing is idempotent", {
  expect_equal(add_unobs_levels(data_full), data_full)
  expect_equal(balance(data_balanced, means), data_balanced)
})

test_that("Re-ordering means doesn't change balancing", {
  shuffled_means <- lapply(means,
                           \(m) m[sample(seq_len(nrow(m))), ])

  expect_equal(balance(data_full, shuffled_means), data_balanced)
})
