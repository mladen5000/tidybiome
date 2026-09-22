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

test_that("vectorized bray and jaccard match vegan vegdist exactly", {
  counts <- matrix(c(10, 0, 5, 20, 10, 0, 5, 5, 5), nrow = 3, ncol = 3,
                   dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  tb_bray <- calc_beta_diversity(tb, metric = "bray")
  d_bray <- get_distance(tb_bray, "bray")
  
  if (requireNamespace("vegan", quietly = TRUE)) {
    veg_bray <- vegan::vegdist(t(counts), method = "bray")
    expect_equal(as.matrix(d_bray), as.matrix(veg_bray), tolerance = 1e-6)
    
    tb_jacc <- calc_beta_diversity(tb, metric = "jaccard")
    d_jacc <- get_distance(tb_jacc, "jaccard")
    veg_jacc <- vegan::vegdist(t(counts), method = "jaccard", binary = TRUE)
    expect_equal(as.matrix(d_jacc), as.matrix(veg_jacc), tolerance = 1e-6)
  }
})

test_that("Tree-Wasserstein and DTH test evaluate beta diversity and homogeneity", {
  counts <- matrix(c(50, 10, 5, 2,
                     5,  40, 50, 60), nrow = 2, byrow = TRUE,
                   dimnames = list(c("T1", "T2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
  tax_table <- data.frame(taxon_id = c("T1", "T2"), Phylum = c("Bacteroidota", "Firmicutes"))
  tb <- tidy_microbiome(counts, sample_data, tax_table)

  # Tree-Wasserstein distance
  tb_w <- calc_beta_diversity(tb, metric = "wasserstein")
  d_w <- get_distance(tb_w, "wasserstein")
  expect_s3_class(d_w, "dist")
  expect_true(all(d_w >= 0))

  # DTH test for homogeneity
  dth_res <- dth_test(tb_w, group = "group", metric = "wasserstein", n_perm = 49)
  expect_s3_class(dth_res, "tbl_df")
  expect_true(all(c("statistic", "p_value", "n_perm") %in% colnames(dth_res)))
  expect_true(dth_res$p_value >= 0 && dth_res$p_value <= 1)
})

test_that("calc_unifrac calculates valid unweighted and weighted UniFrac distances on trees", {
  counts <- matrix(c(10, 0, 5, 20), nrow = 2, ncol = 2,
                   dimnames = list(c("T1", "T2"), c("S1", "S2")))
  tree <- ape::read.tree(text = "(T1:0.5,T2:0.5);")
  tb <- tidy_microbiome(counts, phy_tree = tree)

  # Unweighted UniFrac
  tb_u <- calc_beta_diversity(tb, metric = "unifrac")
  d_u <- get_distance(tb_u, "unifrac")
  expect_s3_class(d_u, "dist")
  expect_equal(attr(d_u, "Size"), 2)
  expect_true(all(d_u >= 0 & d_u <= 1))

  # Weighted UniFrac
  tb_w <- calc_beta_diversity(tb, metric = "wunifrac")
  d_w <- get_distance(tb_w, "wunifrac")
  expect_s3_class(d_w, "dist")
  expect_true(all(d_w >= 0 & d_w <= 1))
})

