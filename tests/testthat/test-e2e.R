library(testthat)
library(tidybiome)

test_that("full end-to-end tidybiome pipeline runs flawlessly", {
  # Load demo data
  data("gut_microbiome", package = "tidybiome")
  expect_s3_class(gut_microbiome, "tidy_microbiome")

  # 1. Tidyverse verbs
  tb_filt <- gut_microbiome |>
    dplyr::filter(depth > 1000) |>
    dplyr::mutate(log_depth = log(depth))
  expect_s3_class(tb_filt, "tidy_microbiome")
  expect_true("log_depth" %in% colnames(tb_filt))

  # 2. Modern transformations (rCLR)
  tb_trans <- transform_abundance(tb_filt, method = "rclr")
  expect_true("rclr" %in% names(attr(tb_trans, "assays")))

  # 3. Modern alpha diversity (Hill numbers)
  tb_alpha <- calc_alpha_diversity(tb_trans, metrics = c("hill", "shannon"))
  expect_true(all(c("hill_0", "hill_1", "hill_2", "shannon") %in% colnames(tb_alpha)))

  # 4. Modern beta diversity (Robust Aitchison)
  tb_beta <- calc_beta_diversity(tb_alpha, metric = "raitchison")
  d <- get_distance(tb_beta, "raitchison")
  expect_s3_class(d, "dist")

  # 5. Modern ordination (RPCA biplot)
  tb_ord <- calc_ordination(tb_beta, method = "rpca")
  ord <- get_ordination(tb_ord, "rpca")
  expect_true(all(c("samples", "taxa", "variance_explained") %in% names(ord)))

  # 6. SOTA Consensus Differential Abundance
  da_res <- calc_differential_abundance(tb_ord, group = "treatment")
  expect_s3_class(da_res, "tbl_df")
  expect_true(all(c("taxon_id", "log2fc", "p_consensus", "padj_consensus", "agreement_score", "is_significant") %in% colnames(da_res)))

  # 7. SOTA Core Microbiome
  core_res <- calc_core_microbiome(tb_ord, method = "threshold", min_prevalence = 0.7)
  expect_s3_class(core_res, "tbl_df")

  # 8. Aesthetic plotting functions
  p1 <- plot_composition(tb_ord, rank = "Phylum", top_n = 4)
  p2 <- plot_ordination(tb_ord, method = "rpca", color = "treatment", biplot = TRUE)
  p3 <- plot_alpha(tb_alpha, metric = "hill_1", x = "treatment")
  p4 <- plot_da_volcano(da_res)
  p5 <- plot_core(core_res)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
})
