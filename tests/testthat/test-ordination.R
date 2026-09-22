library(testthat)
library(tidybiome)

test_that("RPCA generates sample coordinates, taxon loadings, and explained variance", {
  counts <- matrix(c(100, 20, 5, 2,
                     5, 10, 80, 70,
                     1, 1, 50, 40), nrow = 3, byrow = TRUE,
                   dimnames = list(c("Tax1", "Tax2", "Tax3"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
  tax_table <- data.frame(taxon_id = c("Tax1", "Tax2", "Tax3"), Phylum = c("P1", "P2", "P3"))
  tb <- tidy_microbiome(counts, sample_data, tax_table)

  tb_ord <- calc_ordination(tb, method = "rpca")
  ord <- get_ordination(tb_ord, "rpca")

  expect_true(all(c("samples", "taxa", "variance_explained") %in% names(ord)))
  expect_equal(nrow(ord$samples), 4)
  expect_equal(nrow(ord$taxa), 3)
  expect_true(all(c("PC1", "PC2", "group") %in% colnames(ord$samples)))
  expect_true(all(c("PC1", "PC2", "Phylum") %in% colnames(ord$taxa)))
  expect_true(sum(ord$variance_explained) <= 1.0001)
})

test_that("PCoA with negative eigenvalue correction runs smoothly", {
  counts <- matrix(c(10, 0, 5, 20, 10, 0, 5, 5, 5), nrow = 3, ncol = 3,
                   dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  tb_pcoa <- calc_ordination(tb, method = "pcoa", metric = "raitchison")
  ord <- get_ordination(tb_pcoa, "pcoa")

  expect_s3_class(ord$samples, "tbl_df")
  expect_equal(nrow(ord$samples), 3)
  expect_true("PCoA1" %in% colnames(ord$samples))
})

test_that("Contrastive PCA (cPCA) extracts phenotype-specific axes", {
  counts <- matrix(c(100, 110, 5, 2,
                     5, 10, 80, 70,
                     1, 1, 50, 40), nrow = 3, byrow = TRUE,
                   dimnames = list(c("Tax1", "Tax2", "Tax3"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("Control", "Control", "Treated", "Treated"))
  tb <- tidy_microbiome(counts, sample_data)

  tb_cpca <- calc_ordination(tb, method = "cpca", contrast_group = "group", contrast_alpha = 1.0)
  ord <- get_ordination(tb_cpca, "cpca")

  expect_true(all(c("samples", "taxa", "variance_explained") %in% names(ord)))
  expect_true("cPC1" %in% colnames(ord$samples))
  expect_equal(nrow(ord$samples), 4)
})
