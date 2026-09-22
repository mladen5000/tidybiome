library(testthat)
library(tidybiome)

test_that("run_permanova executes vegan::adonis2 and outputs tidy tibble", {
  data("gut_microbiome", package = "tidybiome")
  res <- run_permanova(gut_microbiome, ~ treatment, method = "bray", permutations = 199)
  
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("term", "df", "sum_sq", "r2", "f_stat", "p_value") %in% colnames(res)))
  expect_equal(res$term[1], "treatment")
  expect_true(res$r2[1] > 0 && res$r2[1] <= 1)
  expect_true(res$p_value[1] >= 0 && res$p_value[1] <= 1)
})

test_that("run_betadisper executes vegan::betadisper and returns tidy dispersion metrics", {
  data("gut_microbiome", package = "tidybiome")
  disp <- run_betadisper(gut_microbiome, group = "treatment", method = "bray", permutations = 199)
  
  expect_s3_class(disp, "tidybiome_betadisper")
  expect_s3_class(disp$sample_distances, "tbl_df")
  expect_equal(nrow(disp$sample_distances), nrow(gut_microbiome))
  expect_true(all(c("sample_id", "group", "distance_to_centroid") %in% colnames(disp$sample_distances)))
  expect_s3_class(disp$group_summary, "tbl_df")
  expect_s3_class(disp$test, "tbl_df")
  expect_true(!is.na(disp$test$p_value[1]))
  expect_output(print(disp), "tidybiome_betadisper")
})

test_that("run_nmds computes NMDS coordinates and attaches metadata", {
  data("gut_microbiome", package = "tidybiome")
  nmds_res <- run_nmds(gut_microbiome, assay = "counts", method = "bray", k = 2, trymax = 20)
  
  expect_s3_class(nmds_res, "tbl_df")
  expect_true(all(c("sample_id", "NMDS1", "NMDS2", "treatment") %in% colnames(nmds_res)))
  expect_true(is.numeric(attr(nmds_res, "stress")))
  expect_true(attr(nmds_res, "stress") < 0.3)
})

test_that("run_dbrda performs distance-based redundancy analysis with tidy biplots", {
  data("gut_microbiome", package = "tidybiome")
  dbrda_res <- run_dbrda(gut_microbiome, ~ treatment + age, distance = "bray")
  
  expect_s3_class(dbrda_res, "tidybiome_dbrda")
  expect_s3_class(dbrda_res$samples, "tbl_df")
  expect_true(all(c("sample_id", "dbRDA1", "dbRDA2", "treatment", "age") %in% colnames(dbrda_res$samples)))
  expect_true(is.numeric(dbrda_res$variance_explained))
  expect_true(length(dbrda_res$variance_explained) >= 2)
  expect_output(print(dbrda_res), "tidybiome_dbrda")
})

