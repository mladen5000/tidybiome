test_that("all plotting functions return valid ggplot2 objects", {
  counts <- matrix(c(100, 20, 5, 2,
                     5, 10, 80, 70), nrow = 2, byrow = TRUE,
                   dimnames = list(c("Tax1", "Tax2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"),
                            group = c("A", "A", "B", "B"))
  tax_table <- data.frame(taxon_id = c("Tax1", "Tax2"),
                          Phylum = c("Bacteroidota", "Firmicutes"),
                          Genus = c("Bacteroides", "Faecalibacterium"))
  tb <- tidy_microbiome(counts, sample_data, tax_table)

  # plot_composition
  p_comp <- plot_composition(tb, rank = "Genus", top_n = 1)
  expect_s3_class(p_comp, "ggplot")

  # plot_alpha
  tb <- calc_alpha_diversity(tb, metrics = "hill")
  p_alpha <- plot_alpha(tb, metric = "hill_1", x = "group")
  expect_s3_class(p_alpha, "ggplot")

  # plot_ordination
  tb <- calc_ordination(tb, method = "rpca")
  p_ord <- plot_ordination(tb, method = "rpca", color = "group")
  expect_s3_class(p_ord, "ggplot")

  # plot_core
  core_df <- calc_core_microbiome(tb, method = "threshold", min_prevalence = 0.5)
  p_core <- plot_core(core_df)
  expect_s3_class(p_core, "ggplot")

  # plot_da_volcano
  da_df <- calc_differential_abundance(tb, group = "group")
  p_volc <- plot_da_volcano(da_df)
  expect_s3_class(p_volc, "ggplot")
})
