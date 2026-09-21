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
