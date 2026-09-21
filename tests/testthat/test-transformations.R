test_that("rclr handles zero-inflated matrices without NaNs or pseudocounts", {
  counts <- matrix(c(100, 0, 50, 0,
                     20,  10, 0, 5), nrow = 2, byrow = TRUE,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
  tb <- tidy_microbiome(counts)

  tb_rclr <- transform_abundance(tb, method = "rclr")
  rclr_mat <- attr(tb_rclr, "assays")$rclr
  expect_false(any(is.nan(rclr_mat)))
  expect_false(any(is.infinite(rclr_mat)))

  # Zeros remain zero in robust CLR
  expect_equal(rclr_mat[1, 2], 0)
  expect_equal(rclr_mat[2, 3], 0)

  # Non-zeros are centered by geometric mean of positive values
  pos_vals <- counts[counts[, 1] > 0, 1]
  geom_mean <- exp(mean(log(pos_vals)))
  expect_equal(rclr_mat[1, 1], log(100 / geom_mean))
})

test_that("relabundance scales each sample to sum to 1", {
  counts <- matrix(c(10, 90, 50, 50), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  tb <- tidy_microbiome(counts)
  tb_rel <- transform_abundance(tb, method = "relabundance")
  expect_equal(colSums(attr(tb_rel, "assays")$relabundance), c(S1 = 1, S2 = 1))
})

test_that("clr with pseudocount preserves dimensions and centers each column", {
  counts <- matrix(c(10, 0, 5, 20), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  tb <- tidy_microbiome(counts)
  tb_clr <- transform_abundance(tb, method = "clr", pseudocount = 1)
  clr_mat <- attr(tb_clr, "assays")$clr
  expect_equal(colMeans(clr_mat), c(S1 = 0, S2 = 0))
})
