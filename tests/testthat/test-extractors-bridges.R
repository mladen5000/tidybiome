test_that("tidy_abundance extracts long format and rolls up taxonomy", {
  counts <- matrix(c(10, 5, 20, 8), nrow = 2, ncol = 2,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2")))
  sample_data <- data.frame(sample_id = c("S1", "S2"), group = c("A", "B"))
  tax_table <- data.frame(taxon_id = c("ASV1", "ASV2"),
                          Phylum = c("Bacteroidota", "Bacteroidota"),
                          Genus = c("Bacteroides", "Prevotella"))
  tb <- tidy_microbiome(counts, sample_data, tax_table)

  long_df <- tidy_abundance(tb, assay = "counts", long = TRUE)
  expect_true(all(c("sample_id", "taxon_id", "abundance", "group", "Genus") %in% colnames(long_df)))
  expect_equal(nrow(long_df), 4)

  # Agglomerate by Phylum
  tb_phylum <- aggregate_taxa(tb, rank = "Phylum")
  expect_equal(nrow(attr(tb_phylum, "assays")$counts), 1)
  expect_equal(attr(tb_phylum, "assays")$counts["Bacteroidota", "S1"], 15)
  expect_equal(attr(tb_phylum, "assays")$counts["Bacteroidota", "S2"], 28)
})

test_that("tidy_taxa extracts clean taxonomy table", {
  counts <- matrix(c(10, 5), nrow = 2, ncol = 1,
                   dimnames = list(c("ASV1", "ASV2"), "S1"))
  tax_table <- data.frame(taxon_id = c("ASV1", "ASV2"),
                          Phylum = c("Bacteroidota", "Firmicutes"))
  tb <- tidy_microbiome(counts, tax_table = tax_table)

  tt <- tidy_taxa(tb)
  expect_s3_class(tt, "tbl_df")
  expect_equal(nrow(tt), 2)
  expect_equal(tt$Phylum, c("Bacteroidota", "Firmicutes"))
})
