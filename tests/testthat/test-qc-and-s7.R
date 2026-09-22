test_that("summary.tidy_microbiome produces informative ecological summaries", {
  data(gut_microbiome)
  s <- summary(gut_microbiome)

  expect_s3_class(s, "summary_tidy_microbiome")
  expect_equal(s$n_samples, 24)
  expect_equal(s$n_taxa, 40)
  expect_true(s$sparsity > 0 && s$sparsity < 1)
  expect_true(!is.null(s$depth_summary))
  expect_true(s$depth_summary["Min"] > 0)
  expect_true(s$depth_summary["Max"] >= s$depth_summary["Min"])
  expect_false(s$tree_info$present)
  expect_true("counts" %in% s$assays)
  expect_true("treatment" %in% names(s$metadata_cols))

  # Test print method output without tree
  output <- capture.output(print(s))
  expect_true(any(grepl("tidy_microbiome Summary", output)))
  expect_true(any(grepl("Sequencing Depth", output)))
  expect_true(any(grepl("Taxonomic Hierarchy", output)))
  expect_true(any(grepl("Sample Metadata", output)))

  # Test summary with phylogenetic tree attached
  if (requireNamespace("ape", quietly = TRUE)) {
    taxa_ids <- attr(gut_microbiome, "tax_table")$taxon_id
    tr <- ape::rtree(length(taxa_ids), tip.label = taxa_ids)
    tb_tree <- set_tree(gut_microbiome, tr)
    s_tree <- summary(tb_tree)
    expect_true(s_tree$tree_info$present)
    expect_equal(s_tree$tree_info$n_tips, 40)
    output_tree <- capture.output(print(s_tree))
    expect_true(any(grepl("Phylogenetic Tree", output_tree)))
  }
})

test_that("calc_qc_metrics computes accurate per-sample metrics", {
  data(gut_microbiome)

  # Test augment = TRUE (default)
  tb_qc <- calc_qc_metrics(gut_microbiome, augment = TRUE)
  expect_s3_class(tb_qc, "tidy_microbiome")

  expected_cols <- c("qc_total_reads", "qc_n_features", "qc_sparsity", "qc_top_taxon_share", "qc_shannon")
  expect_true(all(expected_cols %in% colnames(tb_qc)))
  expect_true(all(tb_qc$qc_total_reads > 0))
  expect_true(all(tb_qc$qc_n_features <= 40))
  expect_true(all(tb_qc$qc_sparsity >= 0 & tb_qc$qc_sparsity <= 1))
  expect_true(all(tb_qc$qc_top_taxon_share > 0 & tb_qc$qc_top_taxon_share <= 1))
  expect_true(all(tb_qc$qc_shannon > 0))

  # Test augment = FALSE
  qc_df <- calc_qc_metrics(gut_microbiome, augment = FALSE)
  expect_s3_class(qc_df, "tbl_df")
  expect_equal(nrow(qc_df), 24)
  expect_true("sample_id" %in% colnames(qc_df))
  expect_equal(qc_df$qc_total_reads, tb_qc$qc_total_reads)

  # Test downstream filtering on QC metrics
  tb_filtered <- tb_qc[tb_qc$qc_total_reads >= 5000, ]
  expect_s3_class(tb_filtered, "tidy_microbiome")
  expect_equal(ncol(assay(tb_filtered)), nrow(tb_filtered))
})

test_that("plot_qc creates publication-ready diagnostics", {
  data(gut_microbiome)

  p1 <- plot_qc(gut_microbiome, color_by = "treatment", threshold_x = 5000, threshold_y = 15)
  expect_s3_class(p1, "ggplot")

  p2 <- plot_qc(gut_microbiome, x = "sparsity", y = "top_taxon_share", log_x = FALSE)
  expect_s3_class(p2, "ggplot")
})

test_that("S7 class conversions work bidirectionally when S7 is available", {
  skip_if_not_installed("S7")

  data(gut_microbiome)
  s7_obj <- to_s7(gut_microbiome)

  expect_true(inherits(s7_obj, "tidybiome::tidy_microbiome_s7"))
  meta_s7 <- S7::prop(s7_obj, "metadata")
  expect_equal(nrow(meta_s7), 24)
  expect_true("sample_id" %in% colnames(meta_s7))

  # Test from_s7 roundtrip
  tb_restored <- from_s7(s7_obj)
  expect_s3_class(tb_restored, "tidy_microbiome")
  expect_equal(nrow(tb_restored), nrow(gut_microbiome))
  expect_equal(dim(assay(tb_restored)), dim(assay(gut_microbiome)))
  expect_equal(attr(tb_restored, "tax_table"), attr(gut_microbiome, "tax_table"))

  # Test as_tidybiome generic dispatch
  tb_generic <- as_tidybiome(s7_obj)
  expect_s3_class(tb_generic, "tidy_microbiome")
  expect_equal(nrow(tb_generic), nrow(gut_microbiome))
})
