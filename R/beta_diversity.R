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
                                metric = c("raitchison", "aitchison", "bray", "jaccard", "jsd"),
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
  n <- ncol(mat)
  d_mat <- matrix(0, n, n)
  for (j in 1:(n - 1)) {
    for (k in (j + 1):n) {
      denom <- sum(mat[, j] + mat[, k])
      num <- sum(abs(mat[, j] - mat[, k]))
      val <- if (denom > 0) num / denom else 0
      d_mat[k, j] <- val
      d_mat[j, k] <- val
    }
  }
  dimnames(d_mat) <- list(colnames(mat), colnames(mat))
  stats::as.dist(d_mat)
}

calc_jaccard <- function(mat) {
  pa <- mat > 0
  n <- ncol(mat)
  d_mat <- matrix(0, n, n)
  for (j in 1:(n - 1)) {
    for (k in (j + 1):n) {
      inter <- sum(pa[, j] & pa[, k])
      union_val <- sum(pa[, j] | pa[, k])
      val <- if (union_val > 0) 1 - (inter / union_val) else 0
      d_mat[k, j] <- val
      d_mat[j, k] <- val
    }
  }
  dimnames(d_mat) <- list(colnames(mat), colnames(mat))
  stats::as.dist(d_mat)
}

calc_jsd <- function(mat) {
  # Column-normalize to relative abundance
  prop <- calc_relabundance_matrix(mat)
  n <- ncol(prop)
  d_mat <- matrix(0, n, n)

  for (j in 1:(n - 1)) {
    p <- prop[, j]
    for (k in (j + 1):n) {
      q <- prop[, k]
      m <- 0.5 * (p + q)

      kl_pm <- sum(ifelse(p > 0, p * log(p / m), 0))
      kl_qm <- sum(ifelse(q > 0, q * log(q / m), 0))
      jsd_val <- 0.5 * kl_pm + 0.5 * kl_qm
      val <- sqrt(max(0, jsd_val))
      d_mat[k, j] <- val
      d_mat[j, k] <- val
    }
  }
  dimnames(d_mat) <- list(colnames(mat), colnames(mat))
  stats::as.dist(d_mat)
}
