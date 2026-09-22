library(testthat)
library(tidybiome)

test_that("Hill numbers profile calculates q=0, 1, 2 on effective species scale", {
  # Even community with 4 species: effective species should all be 4
  counts <- matrix(c(25, 25, 25, 25), nrow = 4, ncol = 1,
                   dimnames = list(c("T1", "T2", "T3", "T4"), "S1"))
  tb <- tidy_microbiome(counts)
  tb_alpha <- calc_alpha_diversity(tb, metrics = c("hill", "shannon", "simpson"))

  expect_equal(tb_alpha$hill_0, 4) # Richness
  expect_equal(round(tb_alpha$hill_1, 5), 4) # Exp(Shannon)
  expect_equal(round(tb_alpha$hill_2, 5), 4) # Inverse Simpson
  expect_equal(round(tb_alpha$shannon, 5), round(log(4), 5))
  expect_equal(round(tb_alpha$simpson, 5), 0.75) # 1 - 4 * (1/4)^2 = 0.75
})

test_that("Chao1 handles singletons properly", {
  # 5 observed, 2 singletons, 1 doubleton
  counts <- matrix(c(1, 1, 2, 10, 20), nrow = 5, ncol = 1,
                   dimnames = list(paste0("T", 1:5), "S1"))
  tb <- tidy_microbiome(counts)
  tb_alpha <- calc_alpha_diversity(tb, metrics = c("chao1", "observed"))

  expect_equal(tb_alpha$observed, 5)
  # Chao1 = 5 + (2 * 1) / (2 * 2) = 5.5
  expect_equal(tb_alpha$chao1, 5.5)
})

test_that("Simplex compositional variation calculates Aitchison total variance", {
  # Equal proportions -> variance of log proportions is 0
  counts <- matrix(c(10, 10, 10, 10), nrow = 4, ncol = 1,
                   dimnames = list(paste0("T", 1:4), "S1"))
  tb <- tidy_microbiome(counts)
  tb_alpha <- calc_alpha_diversity(tb, metrics = "simplex_variation")

  expect_true("simplex_variation" %in% colnames(tb_alpha))
  expect_equal(tb_alpha$simplex_variation, 0)

  # Unequal proportions
  counts2 <- matrix(c(100, 10, 1), nrow = 3, ncol = 1,
                    dimnames = list(paste0("T", 1:3), "S1"))
  tb2 <- tidy_microbiome(counts2)
  tb_alpha2 <- calc_alpha_diversity(tb2, metrics = "simplex_variation")
  expect_true(tb_alpha2$simplex_variation > 0)
})
