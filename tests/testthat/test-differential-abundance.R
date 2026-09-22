library(testthat)
library(tidybiome)

test_that("calc_differential_abundance runs consensus engines and flags biomarkers", {
  # Taxon 1 strongly differential between groups A and B
  counts <- matrix(c(100, 120, 110,  5,   8,   6,
                     50,  45,  55,  48,  52,  50), nrow = 2, byrow = TRUE,
                   dimnames = list(c("DiffTaxon", "NullTaxon"), paste0("S", 1:6)))
  sample_data <- data.frame(sample_id = paste0("S", 1:6),
                            group = c("A", "A", "A", "B", "B", "B"))
  tb <- tidy_microbiome(counts, sample_data)

  da_res <- calc_differential_abundance(tb, group = "group", methods = c("consensus", "linda", "clr_linear", "wilcoxon"))
  expect_s3_class(da_res, "tbl_df")
  expect_true(all(c("taxon_id", "log2fc", "p_consensus", "padj_consensus", "agreement_score", "is_significant") %in% colnames(da_res)))

  # DiffTaxon should be detected as significant with high agreement
  diff_row <- da_res[da_res$taxon_id == "DiffTaxon", ]
  expect_true(diff_row$is_significant)
  expect_true(diff_row$agreement_score >= 2)
  expect_true(diff_row$padj_consensus < 0.05)
})

test_that("CAFT engine models zero cells and differential abundance", {
  # Taxon 1 has zero cells in group B, high counts in group A
  counts <- matrix(c(100, 80, 90, 0, 0, 0,
                     20, 25, 22, 21, 24, 20), nrow = 2, byrow = TRUE,
                   dimnames = list(c("ZeroTaxon", "NullTaxon"), paste0("S", 1:6)))
  sample_data <- data.frame(sample_id = paste0("S", 1:6), group = c("A", "A", "A", "B", "B", "B"))
  tb <- tidy_microbiome(counts, sample_data)

  da_caft <- calc_differential_abundance(tb, group = "group", methods = "caft")
  expect_s3_class(da_caft, "tbl_df")
  expect_true(all(c("log2fc_caft", "p_caft", "padj_caft") %in% colnames(da_caft)))

  # ZeroTaxon should have significant p_caft
  z_row <- da_caft[da_caft$taxon_id == "ZeroTaxon", ]
  expect_true(z_row$p_caft < 0.05)
})

test_that("calc_differential_abundance adjusts for covariates using formula", {
  set.seed(42)
  counts <- matrix(rpois(100, lambda = 50), nrow = 10, ncol = 10,
                   dimnames = list(paste0("ASV", 1:10), paste0("S", 1:10)))
  # Inject true signal for ASV1
  counts["ASV1", 1:5] <- counts["ASV1", 1:5] * 5
  
  sample_data <- data.frame(
    sample_id = paste0("S", 1:10),
    treatment = rep(c("A", "B"), each = 5),
    age = rnorm(10, mean = 40, sd = 5)
  )
  tb <- tidy_microbiome(counts, sample_data = sample_data)
  res <- calc_differential_abundance(tb, formula = ~ treatment + age, contrast = "treatment")
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("taxon_id", "log2fc", "p_consensus", "padj_consensus") %in% colnames(res)))
  expect_equal(res$taxon_id[1], "ASV1")
  expect_true(res$padj_consensus[1] < 0.05)
})

test_that("calc_differential_abundance supports continuous target variable", {
  set.seed(42)
  counts <- matrix(rpois(100, lambda = 50), nrow = 10, ncol = 10,
                   dimnames = list(paste0("ASV", 1:10), paste0("S", 1:10)))
  # ASV2 strongly correlates with bmi
  bmi_vals <- seq(20, 35, length.out = 10)
  counts["ASV2", ] <- as.integer(bmi_vals * 8)
  
  sample_data <- data.frame(
    sample_id = paste0("S", 1:10),
    bmi = bmi_vals,
    sex = rep(c("M", "F"), 5)
  )
  tb <- tidy_microbiome(counts, sample_data = sample_data)
  res <- calc_differential_abundance(tb, formula = ~ bmi + sex, contrast = "bmi")
  expect_s3_class(res, "tbl_df")
  expect_true("log2fc" %in% colnames(res))
  expect_true(res$padj_consensus[res$taxon_id == "ASV2"] < 0.05)
})

