p <- 4

idmat <- diag(p)
rownames(idmat) <- paste0("trait", 1:p)
colnames(idmat) <- rownames(idmat)

covmats <- list(
  ind  = idmat,
  sire = 2 * idmat
)

n_ind   <- 5
n_sires <- 2

rtable <- rnesteddata(covmats, c(ind = 5, sire = 2))

test_that("Output table has the right dimensions", {
  expect_equal(dim(rtable), c(n_ind * n_sires, p + 2))
})

test_that("Correct number of distinct levels", {
  expect_equal(length(unique(rtable$ind)), n_ind * n_sires)
  expect_equal(length(unique(rtable$sire)), n_sires)
})
