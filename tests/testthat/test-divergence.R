library(testthat)
library(tidybiome)

test_that("calc_divergence computes dissimilarity against cohort median baseline", {
  data("gut_microbiome", package = "tidybiome")
  
  tb_div <- calc_divergence(gut_microbiome, reference = "median", method = "bray", add_to_metadata = TRUE)
  expect_s3_class(tb_div, "tidy_microbiome")
  expect_true("divergence" %in% colnames(tb_div))
  expect_true(all(tb_div$divergence >= 0))
  expect_equal(length(tb_div$divergence), nrow(gut_microbiome))
})

test_that("calc_divergence computes divergence against a reference group (e.g. Control)", {
  data("gut_microbiome", package = "tidybiome")
  
  tb_div <- calc_divergence(
    gut_microbiome,
    reference = list(treatment = "Control"),
    method = "bray",
    add_to_metadata = TRUE
  )
  expect_true("divergence" %in% colnames(tb_div))
  
  # Samples in Control should on average have lower divergence to Control centroid than Treated
  ctrl_div <- mean(tb_div$divergence[tb_div$treatment == "Control"])
  trt_div <- mean(tb_div$divergence[tb_div$treatment == "Treated"])
  expect_true(ctrl_div < trt_div)
})

test_that("calc_divergence returns a tidy tibble when add_to_metadata is FALSE", {
  data("gut_microbiome", package = "tidybiome")
  
  div_df <- calc_divergence(gut_microbiome, reference = "median", method = "bray", add_to_metadata = FALSE)
  expect_s3_class(div_df, "tbl_df")
  expect_true(all(c("sample_id", "divergence") %in% colnames(div_df)))
  expect_equal(nrow(div_df), nrow(gut_microbiome))
})
