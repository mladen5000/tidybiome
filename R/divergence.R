#' Calculate Ecological Divergence Relative to a Reference State
#'
#' Evaluates the dissimilarity of each sample's community profile against a reference
#' baseline (such as a healthy control group, cohort median profile, or designated baseline sample).
#'
#' @param tb A `tidy_microbiome` object.
#' @param reference Reference baseline to compare against. Can be:
#'   \itemize{
#'     \item `"median"`: Cohort median abundance profile across all samples.
#'     \item `"mean"`: Cohort mean abundance profile.
#'     \item A named list specifying a metadata subset (e.g. `list(treatment = "Control")`).
#'     \item A character string identifying a specific reference `sample_id`.
#'   }
#' @param method Dissimilarity metric: `"bray"` (Bray-Curtis), `"jsd"` (Jensen-Shannon divergence),
#'   or `"euclidean"`. Defaults to `"bray"`.
#' @param assay Assay to evaluate. Defaults to `"counts"`.
#' @param add_to_metadata Logical; if `TRUE`, appends column `divergence` to `tb`'s sample metadata.
#'   If `FALSE`, returns a tidy tibble. Defaults to `TRUE`.
#' @return A `tidy_microbiome` object with `divergence` in metadata (if `add_to_metadata = TRUE`),
#'   or a tibble with `sample_id` and `divergence`.
#' @export
calc_divergence <- function(tb, reference = "median", method = "bray", assay = "counts", add_to_metadata = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in `tb`.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  # Convert to proportions (TSS)
  cs <- colSums(mat)
  cs[cs == 0] <- 1
  prop_mat <- sweep(mat, 2, cs, "/")

  # Construct reference profile
  if (is.character(reference) && length(reference) == 1) {
    if (reference == "median") {
      ref_profile <- apply(prop_mat, 1, stats::median)
    } else if (reference == "mean") {
      ref_profile <- rowMeans(prop_mat)
    } else if (reference %in% colnames(prop_mat)) {
      ref_profile <- prop_mat[, reference]
    } else {
      stop(sprintf("Unknown reference specification '%s'.", reference), call. = FALSE)
    }
  } else if (is.list(reference)) {
    # Subset matching samples from metadata
    meta_df <- tibble::as_tibble(tb)
    match_idx <- rep(TRUE, nrow(meta_df))
    for (nm in names(reference)) {
      if (!nm %in% colnames(meta_df)) {
        stop(sprintf("Variable '%s' from reference list not found in sample metadata.", nm), call. = FALSE)
      }
      match_idx <- match_idx & (meta_df[[nm]] == reference[[nm]])
    }
    if (!any(match_idx)) {
      stop("No samples match the specified reference conditions.", call. = FALSE)
    }
    sub_prop <- prop_mat[, match_idx, drop = FALSE]
    ref_profile <- apply(sub_prop, 1, stats::median)
  } else {
    stop("Invalid `reference` argument.", call. = FALSE)
  }

  # Normalize ref_profile to sum to 1
  sum_ref <- sum(ref_profile)
  if (sum_ref > 0) {
    ref_profile <- ref_profile / sum_ref
  }

  # Compute distance per sample
  n_samples <- ncol(prop_mat)
  div_values <- numeric(n_samples)

  for (j in seq_len(n_samples)) {
    samp <- prop_mat[, j]
    if (method == "bray") {
      div_values[j] <- sum(abs(samp - ref_profile)) / sum(samp + ref_profile)
    } else if (method == "jsd") {
      m <- 0.5 * (samp + ref_profile)
      # Avoid 0 log 0
      nz_s <- samp > 0 & m > 0
      nz_r <- ref_profile > 0 & m > 0
      kl_s <- if (any(nz_s)) sum(samp[nz_s] * log(samp[nz_s] / m[nz_s])) else 0
      kl_r <- if (any(nz_r)) sum(ref_profile[nz_r] * log(ref_profile[nz_r] / m[nz_r])) else 0
      div_values[j] <- sqrt(0.5 * kl_s + 0.5 * kl_r)
    } else if (method == "euclidean") {
      div_values[j] <- sqrt(sum((samp - ref_profile)^2))
    } else {
      stop(sprintf("Unsupported method '%s' in calc_divergence.", method), call. = FALSE)
    }
  }

  if (add_to_metadata) {
    tb$divergence <- div_values
    return(tb)
  }

  tibble::tibble(
    sample_id = colnames(prop_mat),
    divergence = div_values
  )
}
