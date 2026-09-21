test_that("tidy_microbiome initializes and displays correctly", {
  counts <- matrix(c(10, 0, 5, 20, 15, 2, 0, 8), nrow = 2, ncol = 4,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"),
                            group = c("A", "A", "B", "B"))
  tax_table <- data.frame(taxon_id = c("ASV1", "ASV2"),
                          Phylum = c("Bacteroidota", "Firmicutes"),
                          Genus = c("Bacteroides", "Lactobacillus"))
  
  tb <- tidy_microbiome(counts, sample_data, tax_table)
  expect_s3_class(tb, "tidy_microbiome")
  expect_s3_class(tb, "tbl_df")
  expect_equal(nrow(tb), 4)
  expect_equal(attr(tb, "assays")$counts, counts)
  expect_equal(attr(tb, "tax_table")$taxon_id, c("ASV1", "ASV2"))
})

test_that("dplyr verbs synchronize sample metadata and assays", {
  counts <- matrix(c(10, 0, 5, 20, 15, 2, 0, 8), nrow = 2, ncol = 4,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"),
                            group = c("A", "A", "B", "B"),
                            depth = c(15, 22, 15, 10))
  tb <- tidy_microbiome(counts, sample_data)
  
  # Filter
  tb_filt <- dplyr::filter(tb, group == "A")
  expect_equal(nrow(tb_filt), 2)
  expect_equal(colnames(attr(tb_filt, "assays")$counts), c("S1", "S2"))
  
  # Mutate
  tb_mut <- dplyr::mutate(tb, log_depth = log(depth))
  expect_true("log_depth" %in% colnames(tb_mut))
  expect_equal(colnames(attr(tb_mut, "assays")$counts), c("S1", "S2", "S3", "S4"))
  
  # Select
  tb_sel <- dplyr::select(tb, group)
  expect_true(all(c("sample_id", "group") %in% colnames(tb_sel)))
  expect_equal(ncol(attr(tb_sel, "assays")$counts), 4)
  
  # Arrange
  tb_arr <- dplyr::arrange(tb, desc(depth))
  expect_equal(tb_arr$sample_id, c("S2", "S1", "S3", "S4"))
  expect_equal(colnames(attr(tb_arr, "assays")$counts), c("S2", "S1", "S3", "S4"))
  
  # Slice
  tb_slc <- dplyr::slice(tb, 1:2)
  expect_equal(nrow(tb_slc), 2)
  expect_equal(colnames(attr(tb_slc, "assays")$counts), c("S1", "S2"))
})
