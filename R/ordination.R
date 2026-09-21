#' Calculate Dimensionality Reduction and Ordination
#'
#' @description
#' Performs ordination on microbiome data, including modern compositional Robust PCA
#' (RPCA, generating feature loadings for true compositional biplots without distortion)
#' and classical Principal Coordinates Analysis (PCoA with negative eigenvalue correction)
#' or standard PCA.
#'
#' @param tb A `tidy_microbiome` object.
#' @param method Character string specifying the ordination method:
#'   * `"rpca"`: Robust PCA (DEICODE-style SVD on the `rclr` matrix). Computes both
#'     sample coordinates and taxon loading vectors for clean biplots.
#'   * `"pcoa"`: Principal Coordinates Analysis with automatic Lingoes correction for
#'     negative eigenvalues.
#'   * `"pca"`: Standard Principal Component Analysis on centered abundance.
#' @param metric Distance metric to use if `method = "pcoa"`. Defaults to `"raitchison"`.
#' @param assay Character string naming the abundance assay to use. Defaults to `"counts"`.
#' @param n_components Number of ordination axes to retain. Defaults to 3.
#'
#' @return An updated `tidy_microbiome` object with ordination results stored in its metadata.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70, 1, 1, 50, 40), nrow = 3, byrow = TRUE,
#'                  dimnames = list(c("T1", "T2", "T3"), c("S1", "S2", "S3", "S4")))
#' tb <- tidy_microbiome(counts)
#' tb <- calc_ordination(tb, method = "rpca")
#' ord <- get_ordination(tb, "rpca")
calc_ordination <- function(tb,
                            method = c("rpca", "pcoa", "pca"),
                            metric = "raitchison",
                            assay = "counts",
                            n_components = 3) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  method <- match.arg(method)
  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  sample_names <- colnames(mat)
  taxon_names  <- rownames(mat)
  sample_df    <- tibble::as_tibble(tb)
  tax_df       <- attr(tb, "tax_table")

  ord_res <- switch(
    method,
    rpca = run_rpca(mat, sample_df, tax_df, n_components = n_components),
    pcoa = run_pcoa(tb, metric = metric, sample_df = sample_df, n_components = n_components),
    pca  = run_pca(mat, sample_df, tax_df, n_components = n_components)
  )

  attr(tb, "metadata")$ordinations[[method]] <- ord_res
  tb
}

#' Extract Ordination Results from tidy_microbiome
#'
#' @param tb A `tidy_microbiome` object.
#' @param method Optional character string naming the ordination method to extract.
#'   If `NULL`, returns the most recently calculated ordination.
#'
#' @return A list containing:
#'   * `samples`: [tibble::tbl_df] of sample coordinates and sample metadata.
#'   * `taxa`: [tibble::tbl_df] of taxon loadings (if available) and taxonomy.
#'   * `variance_explained`: Numeric vector of variance explained by each axis.
#' @export
get_ordination <- function(tb, method = NULL) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  ord_list <- attr(tb, "metadata")$ordinations
  if (length(ord_list) == 0) {
    stop("No ordinations found in tidy_microbiome. Run `calc_ordination()` first.", call. = FALSE)
  }
  if (is.null(method)) {
    return(ord_list[[length(ord_list)]])
  }
  if (!method %in% names(ord_list)) {
    stop(sprintf("Ordination method '%s' not found. Available: %s",
                 method, paste(names(ord_list), collapse = ", ")), call. = FALSE)
  }
  ord_list[[method]]
}

# Robust PCA implementation
run_rpca <- function(mat, sample_df, tax_df, n_components = 3) {
  # 1. rCLR transform
  rclr_mat <- calc_rclr_matrix(mat)

  # 2. Mean center by feature and sample
  # Transpose so rows are samples, columns are taxa
  X <- t(rclr_mat)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_centered[is.na(X_centered)] <- 0

  # 3. SVD
  k <- min(n_components, nrow(X_centered) - 1, ncol(X_centered) - 1)
  k <- max(1, k)
  svd_fit <- svd(X_centered, nu = k, nv = k)

  # Variance explained
  d_sq <- svd_fit$d^2
  var_exp <- d_sq / sum(d_sq)
  var_exp <- var_exp[seq_len(k)]
  names(var_exp) <- paste0("PC", seq_len(k))

  # Sample coordinates: U * D
  sample_coords <- svd_fit$u %*% diag(svd_fit$d[seq_len(k)], nrow = k, ncol = k)
  colnames(sample_coords) <- paste0("PC", seq_len(k))
  sample_coords_df <- as.data.frame(sample_coords)
  sample_coords_df$sample_id <- rownames(X)

  # Merge with sample metadata
  samples_merged <- merge(sample_coords_df, sample_df, by = "sample_id", sort = FALSE)

  # Taxon loadings: V
  taxa_loadings <- svd_fit$v
  colnames(taxa_loadings) <- paste0("PC", seq_len(k))
  taxa_loadings_df <- as.data.frame(taxa_loadings)
  taxa_loadings_df$taxon_id <- colnames(X)

  # Merge with taxonomy
  if (!is.null(tax_df) && nrow(tax_df) > 0) {
    taxa_merged <- merge(taxa_loadings_df, tax_df, by = "taxon_id", sort = FALSE)
  } else {
    taxa_merged <- taxa_loadings_df
  }

  list(
    samples = tibble::as_tibble(samples_merged),
    taxa = tibble::as_tibble(taxa_merged),
    variance_explained = var_exp
  )
}

# Corrected PCoA implementation
run_pcoa <- function(tb, metric, sample_df, n_components = 3) {
  # Get or compute distance
  dist_list <- attr(tb, "metadata")$distances
  if (!metric %in% names(dist_list)) {
    tb <- calc_beta_diversity(tb, metric = metric)
  }
  d <- get_distance(tb, metric)
  d_mat <- as.matrix(d)
  n <- nrow(d_mat)

  # Gower centering
  H <- diag(n) - (1 / n)
  B <- -0.5 * (H %*% (d_mat^2) %*% H)

  # Eigen decomposition
  eig <- eigen(B, symmetric = TRUE)
  values <- eig$values
  vectors <- eig$vectors

  # Lingoes correction if negative eigenvalues exist
  if (any(values < -1e-8)) {
    min_val <- min(values)
    d_sq_corr <- d_mat^2 + 2 * abs(min_val) * (1 - diag(n))
    B <- -0.5 * (H %*% d_sq_corr %*% H)
    eig <- eigen(B, symmetric = TRUE)
    values <- pmax(0, eig$values)
    vectors <- eig$vectors
  }

  pos_idx <- which(values > 1e-8)
  k <- min(n_components, length(pos_idx))
  k <- max(1, k)

  var_exp <- values[seq_len(k)] / sum(values[pos_idx])
  names(var_exp) <- paste0("PCoA", seq_len(k))

  # Coordinates: V * sqrt(lambda)
  coords <- vectors[, seq_len(k), drop = FALSE] %*% diag(sqrt(values[seq_len(k)]), nrow = k, ncol = k)
  colnames(coords) <- paste0("PCoA", seq_len(k))
  coords_df <- as.data.frame(coords)
  coords_df$sample_id <- rownames(d_mat)

  samples_merged <- merge(coords_df, sample_df, by = "sample_id", sort = FALSE)

  list(
    samples = tibble::as_tibble(samples_merged),
    taxa = NULL,
    variance_explained = var_exp
  )
}

# Standard PCA
run_pca <- function(mat, sample_df, tax_df, n_components = 3) {
  X <- t(mat)
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  k <- min(n_components, ncol(pca_fit$x))

  var_exp <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  var_exp <- var_exp[seq_len(k)]
  names(var_exp) <- paste0("PC", seq_len(k))

  sample_coords <- as.data.frame(pca_fit$x[, seq_len(k), drop = FALSE])
  sample_coords$sample_id <- rownames(sample_coords)
  samples_merged <- merge(sample_coords, sample_df, by = "sample_id", sort = FALSE)

  taxa_loadings <- as.data.frame(pca_fit$rotation[, seq_len(k), drop = FALSE])
  taxa_loadings$taxon_id <- rownames(taxa_loadings)
  taxa_merged <- if (!is.null(tax_df)) merge(taxa_loadings, tax_df, by = "taxon_id", sort = FALSE) else taxa_loadings

  list(
    samples = tibble::as_tibble(samples_merged),
    taxa = tibble::as_tibble(taxa_merged),
    variance_explained = var_exp
  )
}
