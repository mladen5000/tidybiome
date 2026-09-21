library(testthat)
library(tidybiome)

test_that("calc_prevalence calculates prevalence and abundance metrics per taxon", {
  data("gut_microbiome", package = "tidybiome")
  prev_df <- calc_prevalence(gut_microbiome, assay = "counts", detection = 0)
  
  expect_s3_class(prev_df, "tbl_df")
  expect_true(all(c("taxon_id", "prevalence", "prevalence_n", "mean_abundance") %in% colnames(prev_df)))
  expect_equal(nrow(prev_df), 40)
  expect_true(all(prev_df$prevalence >= 0 & prev_df$prevalence <= 1))
  expect_true(all(prev_df$prevalence_n >= 0 & prev_df$prevalence_n <= nrow(gut_microbiome)))
})

test_that("calc_dominant identifies per-sample dominant taxa and adds to metadata", {
  data("gut_microbiome", package = "tidybiome")
  
  # Calculate dominant
  tb_dom <- calc_dominant(gut_microbiome, assay = "counts", add_to_metadata = TRUE)
  expect_s3_class(tb_dom, "tidy_microbiome")
  expect_true("dominant_taxa" %in% colnames(tb_dom))
  expect_equal(length(tb_dom$dominant_taxa), nrow(gut_microbiome))
  
  # Check dominant summary attribute or return when add_to_metadata = FALSE
  dom_summary <- calc_dominant(gut_microbiome, assay = "counts", add_to_metadata = FALSE)
  expect_s3_class(dom_summary, "tbl_df")
  expect_true(all(c("dominant_taxa", "n_samples", "frequency") %in% colnames(dom_summary)))
})

test_that("filter_prevalent filters taxa based on prevalence threshold", {
  data("gut_microbiome", package = "tidybiome")
  
  # Filter to taxa present in at least 70% of samples
  tb_filt <- filter_prevalent(gut_microbiome, min_prevalence = 0.70)
  expect_s3_class(tb_filt, "tidy_microbiome")
  
  n_taxa_filt <- nrow(attr(tb_filt, "tax_table"))
  expect_true(n_taxa_filt < 40)
  expect_true(n_taxa_filt > 0)
  expect_equal(nrow(attr(tb_filt, "assays")$counts), n_taxa_filt)
})
