library(tidybiome)
library(testthat)

test_that("tidy_microbiome rejects invalid count matrices defensively", {
  # Negative counts
  neg_mat <- matrix(c(10, -5, 20, 15), nrow = 2, dimnames = list(c("T1", "T2"), c("S1", "S2")))
  expect_error(tidy_microbiome(neg_mat), "cannot contain negative values")

  # NA values in counts
  na_mat <- matrix(c(10, NA, 20, 15), nrow = 2, dimnames = list(c("T1", "T2"), c("S1", "S2")))
  expect_error(tidy_microbiome(na_mat), "contains NA or NaN values")

  # Inf values in counts
  inf_mat <- matrix(c(10, Inf, 20, 15), nrow = 2, dimnames = list(c("T1", "T2"), c("S1", "S2")))
  expect_error(tidy_microbiome(inf_mat), "cannot contain Infinite values")

  # Empty matrix
  empty_mat <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(tidy_microbiome(empty_mat), "must have at least 1 taxon")

  # Duplicate sample IDs
  dup_s_mat <- matrix(c(10, 5, 20, 15), nrow = 2, dimnames = list(c("T1", "T2"), c("S1", "S1")))
  expect_error(tidy_microbiome(dup_s_mat), "Duplicate column names")

  # Duplicate taxon IDs
  dup_t_mat <- matrix(c(10, 5, 20, 15), nrow = 2, dimnames = list(c("T1", "T1"), c("S1", "S2")))
  expect_error(tidy_microbiome(dup_t_mat), "Duplicate row names")
})

test_that("calc_beta_diversity rejects cohorts with fewer than 2 samples", {
  mat <- matrix(c(10, 20), nrow = 2, ncol = 1, dimnames = list(c("T1", "T2"), "S1"))
  tb <- tidy_microbiome(mat)
  expect_error(calc_beta_diversity(tb, metric = "raitchison"), "requires at least 2 samples")
})

test_that("calc_differential_abundance catches invalid predictors and missing values", {
  counts <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE,
                   dimnames = list(c("T1", "T2"), c("S1", "S2", "S3")))

  # Predictor with single level
  sam_single_lvl <- data.frame(sample_id = c("S1", "S2", "S3"), group = c("A", "A", "A"))
  tb1 <- tidy_microbiome(counts, sample_data = sam_single_lvl)
  expect_error(calc_differential_abundance(tb1, formula = ~ group), "must have at least 2 distinct levels")

  # Continuous predictor with zero variance
  sam_zero_var <- data.frame(sample_id = c("S1", "S2", "S3"), age = c(30, 30, 30))
  tb2 <- tidy_microbiome(counts, sample_data = sam_zero_var)
  expect_error(calc_differential_abundance(tb2, formula = ~ age), "has zero variance")

  # Missing values in formula variables
  sam_na <- data.frame(sample_id = c("S1", "S2", "S3"), group = c("A", "B", NA))
  tb3 <- tidy_microbiome(counts, sample_data = sam_na)
  expect_error(calc_differential_abundance(tb3, formula = ~ group), "contain NA values")
})
