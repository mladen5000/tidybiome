library(testthat)
library(tidybiome)

test_that("set_tree and get_tree manage ape::phylo trees with synchronization", {
  data("gut_microbiome", package = "tidybiome")
  taxa <- attr(gut_microbiome, "tax_table")$taxon_id
  
  # Generate a synthetic tree with ape
  set.seed(42)
  tree <- ape::rtree(n = length(taxa), tip.label = taxa)
  
  # Initially NULL
  expect_null(get_tree(gut_microbiome))
  
  # Attach tree
  tb_with_tree <- set_tree(gut_microbiome, tree)
  expect_s3_class(get_tree(tb_with_tree), "phylo")
  expect_equal(get_tree(tb_with_tree)$tip.label, tree$tip.label)
  
  # Test with tree having extra tips (pruning)
  extra_tree <- ape::rtree(n = length(taxa) + 5, tip.label = c(taxa, paste0("Extra_", 1:5)))
  tb_pruned <- set_tree(gut_microbiome, extra_tree, prune = TRUE)
  expect_equal(length(get_tree(tb_pruned)$tip.label), length(taxa))
})

test_that("as_phyloseq, as_mia, and to_mia aliases are exported", {
  data("gut_microbiome", package = "tidybiome")
  
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    expect_error(as_phyloseq(gut_microbiome), "Package 'phyloseq' is required")
  }
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    expect_error(as_mia(gut_microbiome), "Package 'TreeSummarizedExperiment' is required")
    expect_error(to_mia(gut_microbiome), "Package 'TreeSummarizedExperiment' is required")
  }
})
