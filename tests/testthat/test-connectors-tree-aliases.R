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

test_that("aggregate_taxa preserves tree topology when tree is present", {
  counts <- matrix(c(10, 20, 5, 15), nrow = 2, ncol = 2,
                   dimnames = list(c("T1", "T2"), c("S1", "S2")))
  tax <- data.frame(taxon_id = c("T1", "T2"), Phylum = c("P1", "P2"))
  tree <- ape::read.tree(text = "(T1:1,T2:1);")
  tb <- tidy_microbiome(counts, tax_table = tax, phy_tree = tree)
  
  tb_agg <- aggregate_taxa(tb, rank = "Phylum")
  expect_s3_class(get_tree(tb_agg), "phylo")
  expect_true(all(c("P1", "P2") %in% get_tree(tb_agg)$tip.label))
})

