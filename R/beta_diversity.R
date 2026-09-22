#' Calculate Beta Diversity Dissimilarity and Distance Matrices
#'
#' @description
#' Computes beta diversity distance matrices between samples, including modern
#' compositionally coherent distances (Robust Aitchison, Jensen-Shannon) and classic
#' ecological metrics (Bray-Curtis, Jaccard, classic Aitchison). The resulting
#' distance matrix is stored in the `tidy_microbiome` metadata and can be retrieved
#' using [get_distance()].
#'
#' @param tb A `tidy_microbiome` object.
#' @param metric Character string specifying the distance metric:
#'   * `"raitchison"`: Robust Aitchison distance (Euclidean distance on the `rclr` matrix).
#'     Recommended for zero-inflated microbiome data as it avoids pseudocount distortion.
#'   * `"wasserstein"`: Optimal Transport Tree-Wasserstein / Earth Mover's Distance across
#'     taxonomic hierarchy. Measures physical work needed to move microbial mass across lineages.
#'   * `"aitchison"`: Classic Aitchison distance (Euclidean distance on pseudocount CLR).
#'   * `"bray"`: Bray-Curtis dissimilarity.
#'   * `"jaccard"`: Jaccard binary dissimilarity.
#'   * `"jsd"`: Jensen-Shannon distance (\eqn{\sqrt{\text{JSD}}}).
#' @param assay Character string naming the abundance assay to use. Defaults to `"counts"`.
#' @param pseudocount Numeric pseudocount for `"aitchison"`. Defaults to 1.
#'
#' @return An updated `tidy_microbiome` object with the computed distance stored in its metadata.
#' @export
#'
#' @examples
#' counts <- matrix(c(10, 0, 5, 20, 10, 0, 5, 5, 5), nrow = 3, ncol = 3,
#'                  dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
#' tb <- tidy_microbiome(counts)
#' tb <- calc_beta_diversity(tb, metric = "raitchison")
#' d <- get_distance(tb, "raitchison")
calc_beta_diversity <- function(tb,
                                metric = c("raitchison", "wasserstein", "aitchison", "bray", "jaccard", "jsd"),
                                assay = "counts",
                                pseudocount = 1) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  metric <- match.arg(metric)
  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  sample_names <- colnames(mat)
  n_samp <- length(sample_names)

  dist_obj <- switch(
    metric,
    raitchison = {
      rclr_mat <- calc_rclr_matrix(mat)
      stats::dist(t(rclr_mat), method = "euclidean")
    },
    wasserstein = calc_tree_wasserstein(mat, attr(tb, "tax_table")),
    aitchison = {
      clr_mat <- calc_clr_matrix(mat, pseudocount = pseudocount)
      stats::dist(t(clr_mat), method = "euclidean")
    },
    bray = calc_bray_curtis(mat),
    jaccard = calc_jaccard(mat),
    jsd = calc_jsd(mat)
  )

  attr(dist_obj, "Labels") <- sample_names
  attr(dist_obj, "metric") <- metric

  attr(tb, "metadata")$distances[[metric]] <- dist_obj
  tb
}

#' Extract Distance Matrix from tidy_microbiome
#'
#' @param tb A `tidy_microbiome` object.
#' @param metric Optional character string naming the distance metric to extract.
#'   If `NULL`, returns the most recently calculated distance matrix.
#'
#' @return An object of class `dist`.
#' @export
get_distance <- function(tb, metric = NULL) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  dist_list <- attr(tb, "metadata")$distances
  if (length(dist_list) == 0) {
    stop("No distance matrices found in tidy_microbiome. Run `calc_beta_diversity()` first.", call. = FALSE)
  }
  if (is.null(metric)) {
    return(dist_list[[length(dist_list)]])
  }
  if (!metric %in% names(dist_list)) {
    stop(sprintf("Distance metric '%s' not found. Available: %s",
                 metric, paste(names(dist_list), collapse = ", ")), call. = FALSE)
  }
  dist_list[[metric]]
}

calc_bray_curtis <- function(mat) {
  if (requireNamespace("vegan", quietly = TRUE)) {
    return(vegan::vegdist(t(mat), method = "bray"))
  }
  num <- as.matrix(stats::dist(t(mat), method = "manhattan"))
  cs <- colSums(mat, na.rm = TRUE)
  denom <- outer(cs, cs, "+")
  d_mat <- ifelse(denom > 0, num / denom, 0)
  dimnames(d_mat) <- list(colnames(mat), colnames(mat))
  stats::as.dist(d_mat)
}

calc_jaccard <- function(mat) {
  if (requireNamespace("vegan", quietly = TRUE)) {
    return(vegan::vegdist(t(mat), method = "jaccard", binary = TRUE))
  }
  pa <- (mat > 0) * 1
  inter <- crossprod(pa)
  sums <- colSums(pa)
  union_mat <- outer(sums, sums, "+") - inter
  d_mat <- ifelse(union_mat > 0, 1 - (inter / union_mat), 0)
  diag(d_mat) <- 0
  dimnames(d_mat) <- list(colnames(mat), colnames(mat))
  stats::as.dist(d_mat)
}

calc_jsd <- function(mat) {
  prop <- calc_relabundance_matrix(mat)
  n <- ncol(prop)
  sample_names <- colnames(prop)

  # Precompute Shannon entropy for each column: H(p) = -sum(p * log(p))
  ent <- apply(prop, 2, function(col) {
    pos <- col[col > 0]
    if (length(pos) == 0) return(0)
    -sum(pos * log(pos))
  })

  d_mat <- matrix(0, n, n, dimnames = list(sample_names, sample_names))
  for (j in 1:(n - 1)) {
    p <- prop[, j]
    h_p <- ent[j]
    for (k in (j + 1):n) {
      q <- prop[, k]
      h_q <- ent[k]
      m <- 0.5 * (p + q)
      pos_m <- m[m > 0]
      h_m <- -sum(pos_m * log(pos_m))
      jsd_val <- max(0, h_m - 0.5 * (h_p + h_q))
      val <- sqrt(jsd_val)
      d_mat[k, j] <- val
      d_mat[j, k] <- val
    }
  }
  stats::as.dist(d_mat)
}

calc_tree_wasserstein <- function(mat, tax_df) {
  prop <- calc_relabundance_matrix(mat)
  n_samp <- ncol(prop)
  sample_names <- colnames(prop)

  # Check taxonomy ranks
  rank_cols <- if (!is.null(tax_df)) setdiff(colnames(tax_df), "taxon_id") else character(0)

  # If no ranks, fallback to Total Variation
  if (length(rank_cols) == 0) {
    d_mat <- matrix(0, n_samp, n_samp, dimnames = list(sample_names, sample_names))
    for (j in 1:(n_samp - 1)) {
      for (k in (j + 1):n_samp) {
        val <- 0.5 * sum(abs(prop[, j] - prop[, k]))
        d_mat[k, j] <- val
        d_mat[j, k] <- val
      }
    }
    return(stats::as.dist(d_mat))
  }

  # Build aggregated proportions at each taxonomic level
  level_props <- list()
  weights <- numeric(length(rank_cols) + 1)
  w_val <- 1 / (length(rank_cols) + 1)

  # Leaf level
  level_props[[1]] <- prop
  weights[1] <- w_val

  for (idx in seq_along(rank_cols)) {
    rk <- rank_cols[idx]
    taxa_groups <- as.character(tax_df[[rk]])
    unique_groups <- unique(taxa_groups)
    agg_mat <- matrix(0, nrow = length(unique_groups), ncol = n_samp,
                      dimnames = list(unique_groups, sample_names))
    for (ug in unique_groups) {
      match_idx <- which(taxa_groups == ug)
      if (length(match_idx) == 1) {
        agg_mat[ug, ] <- prop[match_idx, ]
      } else {
        agg_mat[ug, ] <- colSums(prop[match_idx, , drop = FALSE])
      }
    }
    level_props[[idx + 1]] <- agg_mat
    weights[idx + 1] <- w_val
  }

  d_mat <- matrix(0, n_samp, n_samp, dimnames = list(sample_names, sample_names))
  for (j in 1:(n_samp - 1)) {
    for (k in (j + 1):n_samp) {
      w_dist <- 0
      for (l in seq_along(level_props)) {
        mat_l <- level_props[[l]]
        l1_diff <- 0.5 * sum(abs(mat_l[, j] - mat_l[, k]))
        w_dist <- w_dist + weights[l] * l1_diff
      }
      d_mat[k, j] <- w_dist
      d_mat[j, k] <- w_dist
    }
  }
  stats::as.dist(d_mat)
}

#' Distance-based Test for Homogeneity (DTH)
#'
#' @description
#' Implements the 2025 bioRxiv Distance-based Test for Homogeneity (DTH),
#' comparing the empirical distributions of within-group multivariate distances
#' using the 1D Wasserstein optimal transport metric and permutation testing.
#' Unlike `betadisper`, DTH directly evaluates full distance distributions rather
#' than assuming spherical variance around centroids.
#'
#' @param tb A `tidy_microbiome` object.
#' @param group Character string naming the two-level grouping column in sample metadata.
#' @param metric Character string naming the distance metric (default: `"raitchison"`).
#' @param n_perm Number of Monte Carlo permutations (default: 499).
#'
#' @return A [tibble::tbl_df] with `statistic` (Wasserstein distance between dispersion curves),
#'   `p_value`, `n_perm`, and `metric`.
#' @export
dth_test <- function(tb, group, metric = "raitchison", n_perm = 499) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  sample_df <- tibble::as_tibble(tb)
  if (!group %in% colnames(sample_df)) {
    stop(sprintf("Grouping variable '%s' not found.", group), call. = FALSE)
  }

  g_fac <- as.factor(sample_df[[group]])
  lvls <- levels(g_fac)
  if (length(lvls) != 2) {
    stop("`dth_test` currently requires a binary grouping variable with 2 levels.", call. = FALSE)
  }

  dist_list <- attr(tb, "metadata")$distances
  if (!metric %in% names(dist_list)) {
    tb <- calc_beta_diversity(tb, metric = metric)
    dist_list <- attr(tb, "metadata")$distances
  }
  d <- dist_list[[metric]]
  d_mat <- as.matrix(d)

  calc_w1_dispersion <- function(labels) {
    idx1 <- which(labels == lvls[1])
    idx2 <- which(labels == lvls[2])

    sub1 <- d_mat[idx1, idx1][lower.tri(d_mat[idx1, idx1])]
    sub2 <- d_mat[idx2, idx2][lower.tri(d_mat[idx2, idx2])]

    if (length(sub1) == 0 || length(sub2) == 0) return(0)

    # 1D Wasserstein distance between empirical distributions
    probs <- seq(0.01, 0.99, length.out = 50)
    q1 <- stats::quantile(sub1, probs = probs, na.rm = TRUE)
    q2 <- stats::quantile(sub2, probs = probs, na.rm = TRUE)
    mean(abs(q1 - q2))
  }

  obs_stat <- calc_w1_dispersion(g_fac)

  perm_stats <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    shuffled <- sample(g_fac)
    perm_stats[p] <- calc_w1_dispersion(shuffled)
  }

  pval <- (sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)

  tibble::tibble(
    method    = "Distance-based Test for Homogeneity (DTH, bioRxiv 2025)",
    group     = group,
    metric    = metric,
    statistic = obs_stat,
    p_value   = pval,
    n_perm    = n_perm
  )
}
