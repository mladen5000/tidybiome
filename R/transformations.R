#' Transform Abundance Data Using Compositional or Classical Normalization
#'
#' @description
#' Applies modern compositional transformations (such as Robust CLR) or classic
#' normalizations (TSS, pseudocount CLR, presence/absence) to the abundance assay
#' of a `tidy_microbiome` object.
#'
#' @param tb A `tidy_microbiome` object.
#' @param method Character string specifying the transformation method. Options:
#'   * `"rclr"`: Robust Centered Log-Ratio (Martino et al. 2019). Zero values remain zero;
#'     positive values are centered by the geometric mean of observed (non-zero) taxa.
#'   * `"clr"`: Classic Centered Log-Ratio with pseudocount.
#'   * `"relabundance"`: Total Sum Scaling (TSS), scaling each sample to sum to 1.
#'   * `"coverage"`: Coverage-based standardization using Good-Turing sample coverage.
#'   * `"log10"`: Log10 transformation with pseudocount.
#'   * `"pa"`: Presence/Absence binary transformation (0 or 1).
#' @param assay Character string naming the source assay to transform. Defaults to `"counts"`.
#' @param name Optional character string naming the output assay. Defaults to `method`.
#' @param pseudocount Numeric pseudocount added for `"clr"` and `"log10"`. Defaults to 1.
#'
#' @return An updated `tidy_microbiome` object with the new assay added.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 0, 50, 0, 20, 10, 0, 5), nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
#' tb <- tidy_microbiome(counts)
#' tb <- transform_abundance(tb, method = "rclr")
transform_abundance <- function(tb,
                                method = c("rclr", "clr", "relabundance", "coverage", "log10", "pa"),
                                assay = "counts",
                                name = NULL,
                                pseudocount = 1) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  method <- match.arg(method)
  if (is.null(name)) {
    name <- method
  }

  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Source assay '%s' not found in tidy_microbiome.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]

  trans_mat <- switch(
    method,
    rclr = calc_rclr_matrix(mat),
    clr = calc_clr_matrix(mat, pseudocount = pseudocount),
    relabundance = calc_relabundance_matrix(mat),
    coverage = calc_coverage_matrix(mat),
    log10 = calc_log10_matrix(mat, pseudocount = pseudocount),
    pa = calc_pa_matrix(mat)
  )

  attr(tb, "assays")[[name]] <- trans_mat
  tb
}

# Robust CLR: non-zero geometric mean per sample
calc_rclr_matrix <- function(mat) {
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (j in seq_len(ncol(mat))) {
    col_vals <- mat[, j]
    pos_idx <- which(col_vals > 0)
    if (length(pos_idx) > 0) {
      log_pos <- log(col_vals[pos_idx])
      log_geom_mean <- mean(log_pos)
      out[pos_idx, j] <- log_pos - log_geom_mean
    }
  }
  out
}

# Classic CLR with pseudocount
calc_clr_matrix <- function(mat, pseudocount = 1) {
  shifted <- mat + pseudocount
  log_mat <- log(shifted)
  log_geom_means <- colMeans(log_mat)
  sweep(log_mat, 2, log_geom_means, "-")
}

# Relative abundance (TSS)
calc_relabundance_matrix <- function(mat) {
  col_totals <- colSums(mat, na.rm = TRUE)
  col_totals[col_totals == 0] <- 1
  sweep(mat, 2, col_totals, "/")
}

# Coverage-based standardization
calc_coverage_matrix <- function(mat) {
  out <- mat
  for (j in seq_len(ncol(mat))) {
    col_vals <- mat[, j]
    n_reads <- sum(col_vals)
    if (n_reads > 0) {
      f1 <- sum(col_vals == 1)
      coverage <- max(0.01, 1 - (f1 / n_reads))
      out[, j] <- (col_vals / n_reads) * coverage
    }
  }
  out
}

# Log10
calc_log10_matrix <- function(mat, pseudocount = 1) {
  log10(mat + pseudocount)
}

# Presence/Absence
calc_pa_matrix <- function(mat) {
  (mat > 0) * 1
}
