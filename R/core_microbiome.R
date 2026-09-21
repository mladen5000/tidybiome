#' Identify Core Microbiome Members
#'
#' @description
#' Detects the core microbiome members across samples using either modern data-driven
#' inflection curve analysis (identifying natural mathematical breakpoints in the
#' prevalence spectrum) or classic static threshold filtering.
#'
#' @param tb A `tidy_microbiome` object.
#' @param method Character string specifying the core detection method:
#'   * `"inflection"`: Data-driven inflection curve analysis. Evaluates the sorted
#'     prevalence distribution across taxa and finds the natural elbow/inflection point.
#'   * `"threshold"`: Static filter on minimum prevalence and relative abundance.
#' @param min_prevalence Numeric minimum prevalence threshold (between 0 and 1) for `"threshold"` method.
#'   Defaults to 0.8 (present in at least 80% of samples).
#' @param min_abundance Numeric minimum relative abundance threshold for a taxon to be considered present.
#'   Defaults to 0.001 (0.1%).
#' @param assay Character string naming the abundance assay. Defaults to `"counts"`.
#'
#' @return A [tibble::tbl_df] with one row per taxon containing `taxon_id`,
#'   `prevalence`, `mean_abundance`, `is_core`, and taxonomic annotations.
#' @export
#'
#' @examples
#' counts <- matrix(c(50, 40, 60, 55, 45, 0, 0, 10, 0, 0), nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("CoreTaxon", "RareTaxon"), paste0("S", 1:5)))
#' tb <- tidy_microbiome(counts)
#' core_df <- calc_core_microbiome(tb, method = "threshold", min_prevalence = 0.8)
calc_core_microbiome <- function(tb,
                                 method = c("threshold", "inflection"),
                                 min_prevalence = 0.8,
                                 min_abundance = 0.001,
                                 assay = "counts") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  method <- match.arg(method)
  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome.", assay), call. = FALSE)
  }

  counts <- assays[[assay]]
  rel_mat <- calc_relabundance_matrix(counts)
  n_samp <- ncol(rel_mat)

  # Calculate prevalence and mean abundance
  is_present <- rel_mat >= min_abundance
  prev_vec <- rowMeans(is_present)
  mean_abund <- rowMeans(rel_mat)

  taxa_names <- rownames(counts)
  out_df <- tibble::tibble(
    taxon_id       = taxa_names,
    prevalence     = prev_vec,
    mean_abundance = mean_abund
  )

  if (method == "threshold") {
    out_df$is_core <- out_df$prevalence >= min_prevalence
  } else {
    # Data-driven inflection point
    sorted_prev <- sort(prev_vec, decreasing = TRUE)
    n_taxa <- length(sorted_prev)

    if (n_taxa <= 3) {
      out_df$is_core <- out_df$prevalence >= min_prevalence
    } else {
      # First and second differences
      diff1 <- diff(sorted_prev)
      # Find maximum drop in prevalence
      inflection_idx <- which.min(diff1)
      cutoff_val <- sorted_prev[inflection_idx]
      out_df$is_core <- out_df$prevalence >= cutoff_val
    }
  }

  # Merge taxonomy
  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df) && nrow(tax_df) > 0) {
    out_df <- merge(out_df, tax_df, by = "taxon_id", all.x = TRUE, sort = FALSE)
  }

  out_df <- out_df[order(-out_df$prevalence, -out_df$mean_abundance), ]
  out_df <- tibble::as_tibble(out_df)

  attr(tb, "metadata")$core_microbiome <- out_df
  out_df
}
