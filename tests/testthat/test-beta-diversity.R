test_that("Robust Aitchison distance computes Euclidean distance on rclr", {
  counts <- matrix(c(10, 0, 5,
                     20, 10, 0,
                     5,  5,  5), nrow = 3, ncol = 3,
                   dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  tb_dist <- calc_beta_diversity(tb, metric = "raitchison")

  d <- get_distance(tb_dist, "raitchison")
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 3)
  expect_true(all(d >= 0))
  expect_equal(attr(d, "Labels"), c("S1", "S2", "S3"))
})

test_that("Bray-Curtis and JSD distances calculate valid dissimilarities", {
  counts <- matrix(c(10, 20, 0, 5,
                     15, 25, 5, 0), nrow = 2, byrow = TRUE,
                   dimnames = list(c("A", "B"), c("S1", "S2", "S3", "S4")))
  tb <- tidy_microbiome(counts)

  tb_bray <- calc_beta_diversity(tb, metric = "bray")
  d_bray <- get_distance(tb_bray, "bray")
  expect_true(all(d_bray >= 0 & d_bray <= 1))

  tb_jsd <- calc_beta_diversity(tb, metric = "jsd")
  d_jsd <- get_distance(tb_jsd, "jsd")
  expect_true(all(d_jsd >= 0))
})
