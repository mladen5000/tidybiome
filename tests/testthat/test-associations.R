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

test_that("calc_network computes microbial co-occurrence networks with tidy edges and nodes", {
  data("gut_microbiome", package = "tidybiome")
  net <- calc_network(gut_microbiome, method = "spearman", min_prevalence = 0.3, r_cutoff = 0.3, p_cutoff = 0.1)
  
  expect_s3_class(net, "tidybiome_network")
  expect_s3_class(net$nodes, "tbl_df")
  expect_s3_class(net$edges, "tbl_df")
  expect_true(all(c("taxon_id", "degree", "mean_abundance") %in% colnames(net$nodes)))
  if (nrow(net$edges) > 0) {
    expect_true(all(c("from", "to", "correlation", "p_value", "padj", "weight", "direction") %in% colnames(net$edges)))
    expect_true(all(net$edges$weight >= 0.3))
  }
  expect_output(print(net), "tidybiome_network")
})

