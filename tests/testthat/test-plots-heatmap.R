library(tidybiome)
library(testthat)
library(ggplot2)

test_that("plot_heatmap generates valid ggplot objects with clustering and scaling", {
  data("gut_microbiome", package = "tidybiome")
  tb <- gut_microbiome

  p1 <- plot_heatmap(tb, top_n = 10, scale = "log10")
  expect_s3_class(p1, "ggplot")

  p2 <- plot_heatmap(tb, rank = "Genus", top_n = 8, scale = "relabundance", annotation_col = "treatment")
  expect_s3_class(p2, "ggplot")

  p3 <- plot_heatmap(tb, top_n = 5, scale = "rclr", cluster_samples = FALSE, cluster_taxa = FALSE)
  expect_s3_class(p3, "ggplot")
})
