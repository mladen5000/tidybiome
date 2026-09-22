library(tidybiome)
library(testthat)

test_that("GMPR size factor calculation works as expected", {
  mat <- matrix(c(100, 200, 0, 50,
                  50,  100, 10, 25,
                  10,  20,  0, 5,
                  0,   0,   50, 0), nrow = 4, byrow = TRUE,
                dimnames = list(c("T1", "T2", "T3", "T4"), c("S1", "S2", "S3", "S4")))

  sf <- calc_gmpr_size_factors(mat, min_overlap = 2)
  expect_type(sf, "double")
  expect_equal(length(sf), 4)
  expect_true(all(sf > 0))
  # S2 has twice the counts of S1 for overlapping taxa, so sf[S2] should be approx 2 * sf[S1]
  expect_equal(sf[["S2"]] / sf[["S1"]], 2, tolerance = 1e-4)
})

test_that("transform_abundance supports gmpr and hellinger", {
  counts <- matrix(c(100, 200, 50,
                     50,  100, 25,
                     10,  20,  5), nrow = 3, byrow = TRUE,
                   dimnames = list(c("T1", "T2", "T3"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)

  tb <- transform_abundance(tb, method = "gmpr")
  expect_true("gmpr" %in% names(attr(tb, "assays")))
  expect_equal(dim(assay(tb, "gmpr")), c(3, 3))

  tb <- transform_abundance(tb, method = "hellinger")
  expect_true("hellinger" %in% names(attr(tb, "assays")))
  h_mat <- assay(tb, "hellinger")
  # For each sample in hellinger, sum of squares should be 1
  expect_equal(as.numeric(colSums(h_mat^2)), rep(1, 3), tolerance = 1e-6)
})
