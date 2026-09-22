#' Calculate Alpha Diversity and Simplex Compositional Metrics
#'
#' @description
#' Calculates alpha diversity indices including the unified Hill numbers profile
#' (\eqn{q = 0, 1, 2}), simplex compositional variation (bioRxiv 2026),
#' and classic ecological indices (Shannon, Gini-Simpson, Inverse Simpson,
#' Observed, Chao1, and Faith's PD).
#'
#' @param tb A `tidy_microbiome` object.
#' @param metrics Character vector of metrics to calculate. Supported options:
#'   * `"hill"`: Calculates the Hill numbers profile: `hill_0` (species richness),
#'     `hill_1` (exponential of Shannon entropy), and `hill_2` (inverse Simpson index).
#'   * `"simplex_variation"`: Compositional total variance on the closed simplex (bioRxiv 2026),
#'     measuring true geometric dispersion among taxon ratios.
#'   * `"shannon"`: Shannon diversity index (\eqn{H = -\sum p_i \ln p_i}).
#'   * `"simpson"`: Gini-Simpson index (\eqn{1 - \sum p_i^2}).
#'   * `"inv_simpson"`: Inverse Simpson index (\eqn{1 / \sum p_i^2}).
#'   * `"observed"`: Number of observed taxa with counts > 0.
#'   * `"chao1"`: Chao1 bias-corrected richness estimator.
#'   * `"faith_pd"`: Faith's phylogenetic diversity (requires `phy_tree`).
#' @param assay Character string naming the abundance assay to use. Defaults to `"counts"`.
#'
#' @return An updated `tidy_microbiome` object with new alpha diversity columns added.
#' @export
calc_alpha_diversity <- function(tb,
                                 metrics = c("hill", "simplex_variation", "shannon", "simpson", "inv_simpson", "observed", "chao1"),
                                 assay = "counts") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  n_samp <- ncol(mat)
  sample_names <- colnames(mat)

  calc_single_sample <- function(vec) {
    pos_idx <- which(vec > 0)
    pos_vals <- vec[pos_idx]
    s_obs <- length(pos_vals)
    total <- sum(pos_vals)

    if (total == 0 || s_obs == 0) {
      return(list(
        hill_0 = 0, hill_1 = 0, hill_2 = 0,
        shannon = 0, simpson = 0, inv_simpson = 0,
        observed = 0, chao1 = 0
      ))
    }

    p <- pos_vals / total
    h <- -sum(p * log(p))
    d <- sum(p^2)

    # Chao1
    f1 <- sum(pos_vals == 1)
    f2 <- sum(pos_vals == 2)
    chao1_val <- if (f2 > 0) {
      s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    } else {
      s_obs + (f1 * (f1 - 1)) / 2
    }

    # Simplex compositional variation (Aitchison total variance, bioRxiv 2026)
    simplex_var <- if (s_obs >= 2) {
      log_p <- log(pos_vals)
      stats::var(log_p) * ((s_obs - 1) / s_obs)
    } else {
      0
    }

    list(
      hill_0 = s_obs,
      hill_1 = exp(h),
      hill_2 = if (d > 0) 1 / d else 0,
      simplex_variation = simplex_var,
      shannon = h,
      simpson = 1 - d,
      inv_simpson = if (d > 0) 1 / d else 0,
      observed = s_obs,
      chao1 = chao1_val
    )
  }

  res_list <- lapply(seq_len(n_samp), function(j) calc_single_sample(mat[, j]))

  # Append selected metrics
  sample_df <- tibble::as_tibble(tb)

  for (m in metrics) {
    if (m == "hill") {
      sample_df$hill_0 <- vapply(res_list, function(r) r$hill_0, numeric(1))
      sample_df$hill_1 <- vapply(res_list, function(r) r$hill_1, numeric(1))
      sample_df$hill_2 <- vapply(res_list, function(r) r$hill_2, numeric(1))
    } else if (m == "faith_pd") {
      ptree <- attr(tb, "phy_tree")
      if (is.null(ptree)) {
        warning("`phy_tree` not found in tidy_microbiome. Skipping `faith_pd`.", call. = FALSE)
      } else {
        sample_df$faith_pd <- calc_faith_pd(mat, ptree)
      }
    } else if (m %in% names(res_list[[1]])) {
      sample_df[[m]] <- vapply(res_list, function(r) r[[m]], numeric(1))
    }
  }

  reconstruct_tidy_microbiome(sample_df, tb, sync_samples = FALSE)
}

calc_faith_pd <- function(mat, tree) {
  n_samp <- ncol(mat)
  common_taxa <- intersect(rownames(mat), tree$tip.label)

  if (length(common_taxa) <= 1 || !requireNamespace("ape", quietly = TRUE)) {
    return(numeric(n_samp))
  }

  tree_sub <- ape::keep.tip(tree, common_taxa)
  mat_sub <- mat[tree_sub$tip.label, , drop = FALSE]

  n_edge <- nrow(tree_sub$edge)
  n_tip <- length(tree_sub$tip.label)
  edge_lengths <- tree_sub$edge.length
  edge_lengths[is.na(edge_lengths)] <- 0
  parent <- tree_sub$edge[, 1]
  child  <- tree_sub$edge[, 2]

  tip_to_edge <- vector("list", n_tip)
  for (e_idx in seq_len(n_edge)) {
    c_node <- child[e_idx]
    if (c_node <= n_tip) tip_to_edge[[c_node]] <- c(tip_to_edge[[c_node]], e_idx)
  }
  for (t_idx in seq_len(n_tip)) {
    curr <- child[tip_to_edge[[t_idx]][1]]
    repeat {
      parent_edge <- which(child == parent[which(child == curr)[1]])
      if (length(parent_edge) > 0) {
        tip_to_edge[[t_idx]] <- c(tip_to_edge[[t_idx]], parent_edge)
        curr <- child[parent_edge]
      } else break
    }
  }

  M <- matrix(0, nrow = n_edge, ncol = n_tip)
  for (t_idx in seq_len(n_tip)) M[tip_to_edge[[t_idx]], t_idx] <- 1

  # Active edges across all samples simultaneously
  edge_occ_counts <- M %*% (mat_sub > 0)
  pd_vals <- colSums(edge_lengths * (edge_occ_counts > 0))

  # Samples with <= 1 taxon have 0 PD
  n_present <- colSums(mat_sub > 0)
  pd_vals[n_present <= 1] <- 0

  as.numeric(pd_vals)
}
