library(testthat)
library(tidybiome)

test_that("calc_cross_association computes taxon-metadata correlations with FDR adjustment", {
  data("gut_microbiome", package = "tidybiome")
  
  # Calculate cross-associations for numeric metadata columns (age, depth)
  assoc_df <- calc_cross_association(
    gut_microbiome,
    variables = c("age", "depth"),
    method = "spearman",
    p_adj_method = "BH"
  )
  
  expect_s3_class(assoc_df, "tbl_df")
  expect_true(all(c("taxon_id", "variable", "correlation", "p_value", "padj") %in% colnames(assoc_df)))
  expect_equal(nrow(assoc_df), 40 * 2) # 40 taxa x 2 variables
  expect_true(all(assoc_df$correlation >= -1 & assoc_df$correlation <= 1))
  expect_true(all(assoc_df$p_value >= 0 & assoc_df$p_value <= 1))
  expect_true(all(assoc_df$padj >= 0 & assoc_df$padj <= 1))
})

test_that("calc_cross_association automatically detects numeric variables when none specified", {
  data("gut_microbiome", package = "tidybiome")
  
  assoc_auto <- calc_cross_association(gut_microbiome, method = "pearson")
  expect_s3_class(assoc_auto, "tbl_df")
  # gut_microbiome has numeric columns 'age' and 'depth'
  expect_true(all(c("age", "depth") %in% unique(assoc_auto$variable)))
})
