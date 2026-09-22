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

test_that("bracket subsetting slices sample metadata and assay matrices synchronously", {
  counts <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("T1", "T2"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  
  # Row slicing
  sub_tb <- tb[1:2, ]
  expect_s3_class(sub_tb, "tidy_microbiome")
  expect_equal(nrow(sub_tb), 2)
  expect_equal(colnames(attr(sub_tb, "assays")$counts), c("S1", "S2"))
  expect_equal(sub_tb$sample_id, c("S1", "S2"))
  
  # Column slicing
  sub_col <- tb[, "depth"]
  expect_s3_class(sub_col, "tidy_microbiome")
  expect_true(all(c("sample_id", "depth") %in% colnames(sub_col)))
  expect_equal(ncol(attr(sub_col, "assays")$counts), 3)
})

test_that("filter_taxa slices features across assays, taxonomy, and tree", {
  counts <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("T1", "T2"), c("S1", "S2", "S3")))
  tax <- tibble::tibble(taxon_id = c("T1", "T2"), Phylum = c("Bacteroidota", "Firmicutes"))
  tree <- ape::read.tree(text = "(T1:1,T2:1);")
  tb <- tidy_microbiome(counts, tax_table = tax, phy_tree = tree)
  
  filt_tb <- filter_taxa(tb, Phylum == "Bacteroidota")
  expect_s3_class(filt_tb, "tidy_microbiome")
  expect_equal(nrow(attr(filt_tb, "assays")$counts), 1)
  expect_equal(rownames(attr(filt_tb, "assays")$counts), "T1")
  expect_equal(attr(filt_tb, "tax_table")$taxon_id, "T1")
  expect_equal(attr(filt_tb, "phy_tree")$tip.label, "T1")
})

test_that("validate_tidy_microbiome verifies data integrity", {
  counts <- matrix(1:4, nrow = 2, ncol = 2, dimnames = list(c("T1", "T2"), c("S1", "S2")))
  tb <- tidy_microbiome(counts)
  expect_true(validate_tidy_microbiome(tb))
  
  # Corrupt sample_id order in assays
  bad_tb <- tb
  attr(bad_tb, "assays")$counts <- counts[, c(2, 1)]
  expect_error(validate_tidy_microbiome(bad_tb), "Sample names in assay")
})

