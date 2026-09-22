library(testthat)
library(tidybiome)

test_that("calc_core_microbiome detects core members accurately", {
  counts <- matrix(c(50, 40, 60, 55, 45,
                     0,  0,  10, 0,  0), nrow = 2, byrow = TRUE,
                   dimnames = list(c("CoreTaxon", "RareTaxon"), paste0("S", 1:5)))
  tb <- tidy_microbiome(counts)

  core_res <- calc_core_microbiome(tb, method = "threshold", min_prevalence = 0.8)
  expect_s3_class(core_res, "tbl_df")
  expect_true(all(c("taxon_id", "prevalence", "mean_abundance", "is_core") %in% colnames(core_res)))
  expect_true(core_res$is_core[core_res$taxon_id == "CoreTaxon"])
  expect_false(core_res$is_core[core_res$taxon_id == "RareTaxon"])
})
