#' Calculate Alpha Diversity Indices
#'
#' @description
#' Calculates alpha diversity indices including the unified Hill numbers profile
#' (\eqn{q = 0, 1, 2}) and classic ecological indices (Shannon, Gini-Simpson, Inverse Simpson,
#' Observed, Chao1, and Faith's PD). Calculated values are added directly as columns
#' to the sample metadata of the `tidy_microbiome` object.
#'
#' @param tb A `tidy_microbiome` object.
#' @param metrics Character vector of metrics to calculate. Supported options:
#'   * `"hill"`: Calculates the Hill numbers profile: `hill_0` (species richness),
#'     `hill_1` (exponential of Shannon entropy), and `hill_2` (inverse Simpson index).
#'     All three share the identical intuitive unit: "effective number of species".
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
#'
#' @examples
#' counts <- matrix(c(25, 25, 25, 25), nrow = 4, ncol = 1,
#'                  dimnames = list(c("T1", "T2", "T3", "T4"), "S1"))
#' tb <- tidy_microbiome(counts)
#' tb <- calc_alpha_diversity(tb, metrics = "hill")
calc_alpha_diversity <- function(tb,
                                 metrics = c("hill", "shannon", "simpson", "inv_simpson", "observed", "chao1"),
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

    list(
      hill_0 = s_obs,
      hill_1 = exp(h),
      hill_2 = if (d > 0) 1 / d else 0,
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
  # Simple vectorized Faith's PD
  n_samp <- ncol(mat)
  pd_vals <- numeric(n_samp)
  taxa_names <- rownames(mat)

  for (j in seq_len(n_samp)) {
    present_taxa <- taxa_names[mat[, j] > 0]
    if (length(present_taxa) <= 1) {
      pd_vals[j] <- 0
    } else {
      # Keep tips present in tree
      common_tips <- intersect(present_taxa, tree$tip.label)
      if (length(common_tips) <= 1) {
        pd_vals[j] <- 0
      } else if (requireNamespace("ape", quietly = TRUE)) {
        sub_tree <- ape::keep.tip(tree, common_tips)
        pd_vals[j] <- sum(sub_tree$edge.length, na.rm = TRUE)
      } else {
        pd_vals[j] <- length(common_tips)
      }
    }
  }
  pd_vals
}
